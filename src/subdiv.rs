use crate::{
    iterator::HasIterators,
    topol::{HasTopology, Topology},
    EditableTopology, Error, FloatScalarAdaptor, Handle, PolyMeshT, FH, HH,
};
use core::f64;
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul},
};

fn check_for_deleted(topol: &Topology) -> Result<(), Error> {
    // Ensure no deleted topology.
    let fstatus = topol.fstatus.try_borrow()?;
    if let Some(f) = topol
        .faces()
        .find(|f| fstatus[f.index() as usize].deleted())
    {
        return Err(Error::DeletedFace(f));
    }
    let estatus = topol.estatus.try_borrow()?;
    if let Some(e) = topol
        .edges()
        .find(|e| estatus[e.index() as usize].deleted())
    {
        return Err(Error::DeletedEdge(e));
    }
    let vstatus = topol.hstatus.try_borrow()?;
    if let Some(v) = topol
        .vertices()
        .find(|v| vstatus[v.index() as usize].deleted())
    {
        return Err(Error::DeletedVertex(v));
    }
    Ok(())
}

/// This struct doesn't conatin any data. This purpose of this struct is to
/// provide a convenient scope inside its `impl` where I can impose all the
/// trait bounds once and write several free functions within that scope, with
/// those trait bounds applied.
struct CatmullClark<const DIM: usize, A>(PhantomData<A>);

impl<const DIM: usize, A> CatmullClark<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    /// Reserve the space for the final mesh after `niter` iterations of
    /// subdivision. Also reserve space for buffers used to store the points
    /// corresponding to vertices, edges and faces. These buffers are used to
    /// store the points corresponding to the mesh before subdivision, so the
    /// capacity of these buffers will correspond to one less iteration of
    /// subdivision (hence less memory required) than that of the mesh
    /// itself. This reserve takes care of that.
    fn reserve(
        niter: usize,
        mesh: &mut Topology,
        vertex_points: &mut Option<Vec<A::Vector>>,
        edge_points: &mut Vec<A::Vector>,
        face_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        debug_assert!(niter > 0);
        let mut nv = mesh.num_vertices();
        let mut ne = mesh.num_edges();
        let mut nf = mesh.num_faces();
        for i in 0..niter {
            if i == niter - 1 {
                // Reserve the pos arrays.
                if let Some(vpos) = vertex_points {
                    vpos.reserve(nv);
                }
                edge_points.reserve(ne);
                face_points.reserve(nf);
            }
            let v = nv + ne + nf;
            let f = if i == 0 {
                mesh.faces().map(|f| f.valence(mesh)).sum::<usize>()
            } else {
                nf * 4
            };
            let e = 2 * ne + f;
            (nv, ne, nf) = (v, e, f);
        }
        mesh.reserve(nv, ne, nf)
    }

    /// Compute the locations of the points to split the faces and edges.
    fn calc_face_edge_points(
        mesh: &PolyMeshT<DIM, A>,
        update_points: bool,
        face_points: &mut Vec<A::Vector>,
        edge_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        // Compute face points.
        face_points.clear();
        face_points.extend(mesh.faces().map(|f| mesh.calc_face_centroid(f, &points)));
        // Compute edge points.
        edge_points.clear();
        edge_points.extend(mesh.edges().map(|e| {
            let (h, oh) = e.halfedges();
            match (h.face(mesh), oh.face(mesh)) {
                (Some(fa), Some(fb)) if update_points => {
                    (points[h.head(mesh).index() as usize]
                        + points[oh.head(mesh).index() as usize]
                        + face_points[fa.index() as usize]
                        + face_points[fb.index() as usize])
                        * A::scalarf64(0.25)
                }
                _ => {
                    (points[h.head(mesh).index() as usize] + points[oh.head(mesh).index() as usize])
                        * A::scalarf64(0.5)
                }
            }
        }));
        Ok(())
    }

    /// Compute the new locations of the vertices, and move them there. This
    /// function makes use of the preallocated buffer `vertex_points` to avoid
    /// repeated allocations.
    fn update_vertex_positions(
        mesh: &mut PolyMeshT<DIM, A>,
        face_points: &[A::Vector],
        edge_points: &[A::Vector],
        vertex_points: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        // Compute vertex positions.
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut()?;
        vertex_points.clear();
        {
            let points: &[A::Vector] = &points;
            vertex_points.extend(mesh.vertices().map(|v| {
                if v.is_boundary(mesh) {
                    let (count, sum) = mesh
                        .ve_ccw_iter(v)
                        .filter(|e| e.is_boundary(mesh))
                        .fold((1usize, points[v.index() as usize]), |(count, total), e| {
                            (count + 1, total + edge_points[e.index() as usize])
                        });
                    sum / A::scalarf64(count as f64)
                } else {
                    let valence = v.valence(mesh) as f64;
                    (((mesh.vf_ccw_iter(v).fold(A::zero_vector(), |total, f| {
                        total + face_points[f.index() as usize]
                    }) + mesh.vv_ccw_iter(v).fold(A::zero_vector(), |total, v| {
                        total + points[v.index() as usize]
                    })) / A::scalarf64(valence))
                        + (points[v.index() as usize] * A::scalarf64(valence - 2.0)))
                        / A::scalarf64(valence)
                }
            }));
        }
        // Update the vertex positions.
        points.copy_from_slice(vertex_points);
        Ok(())
    }

    /// Split a face according to the catmull clark scheme. This must be called
    /// after all the edges of the face are already split. `hloop`, `spliths`,
    /// and `subfaces` are preallocated buffers used by this function to avoid
    /// repeated allocations inside the call.
    fn split_face(
        mesh: &mut PolyMeshT<DIM, A>,
        f: FH,
        num_old_verts: u32,
        face_points: &[A::Vector],
        hloop: &mut Vec<HH>,
        spliths: &mut Vec<HH>,
        subfaces: &mut Vec<FH>,
    ) -> Result<(), Error> {
        // Find a halfedge that points to an old vertex, and collect the
        // loop of halfedges starting from there.
        let hstart = mesh
            .fh_ccw_iter(f)
            .find(|&h| h.head(mesh).index() < num_old_verts)
            .ok_or(Error::CannotSplitFace(f))?;
        hloop.clear();
        hloop.extend(mesh.loop_ccw_iter(hstart));
        let hloop: &[HH] = hloop; // Immutable.
        debug_assert!(hloop.len() % 2 == 0);
        let valence = mesh.loop_ccw_iter(hstart).count() / 2;
        debug_assert_eq!(valence * 2, hloop.len());
        let ne = mesh.num_edges();
        // New vertex in the middle.
        let fv = mesh.add_vertex(face_points[f.index() as usize])?;
        // Create new edges and faces, with some math to get the indices
        // of edges we haven't added yet. This is a bit sketchy, but
        // should be safe with ample testing.
        spliths.clear();
        subfaces.clear();
        for (lei, hpair) in hloop.chunks_exact(2).enumerate() {
            let h1 = hpair[0];
            let pei = (ne + ((lei + valence - 1) % valence)) as u32;
            let nei = (ne + ((lei + 1) % valence)) as u32;
            let enew = mesh.topol.new_edge(
                h1.tail(mesh),
                fv,
                h1.prev(mesh),
                (2 * pei + 1).into(),
                (2 * nei).into(),
                h1,
            )?;
            debug_assert_eq!(enew.index(), (lei + ne) as u32);
            spliths.push(enew.halfedge(false));
            let flocal = if h1 == hstart {
                f
            } else {
                mesh.topol.new_face(h1)?
            };
            subfaces.push(flocal);
        }
        for (i, hpair) in hloop.chunks_exact(2).enumerate() {
            let rh = spliths[i];
            let orh = rh.opposite();
            let flocal = subfaces[i];
            let pflocal = subfaces[(i + valence - 1) % valence];
            let h1 = hpair[0];
            let h2 = hpair[1];
            // Link halfedges and faces.
            mesh.topol.face_mut(flocal).halfedge = h1;
            mesh.topol.halfedge_mut(rh).face = Some(pflocal);
            mesh.topol.halfedge_mut(orh).face = Some(flocal);
            mesh.topol.halfedge_mut(h1).face = Some(flocal);
            mesh.topol.halfedge_mut(h2).face = Some(flocal);
            // Link halfedges.
            mesh.topol.link_halfedges(rh.prev(mesh), rh);
            mesh.topol.link_halfedges(orh, h1);
            mesh.topol.link_halfedges(h1, h2);
        }
        mesh.topol.vertex_mut(fv).halfedge = Some(spliths[0].opposite());
        Ok(())
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    /// Subdivide the mesh according to the [Catmull-Clark
    /// scheme](https://en.wikipedia.org/wiki/Catmull%E2%80%93Clark_subdivision_surface).
    ///
    /// Subdivisions are carried out for the given number of
    /// `iterations`. `update_points` determines whether the vertices of the
    /// mesh are moved. If this is `false`, only the topology of the mesh is
    /// subdivided, leaving the shape of the unmodified. If this is `true`, then
    /// the shape is smoothed according to the Catmull-Clark scheme.
    ///
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
    /// assert_eq!((8, 12, 6), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.subdivide_catmull_clark(1, true)
    ///     .expect("Subdivision failed");
    /// // The mesh now has more faces.
    /// assert_eq!((26, 48, 24), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.check_topology().expect("Topological errors found");
    /// ```
    pub fn subdivide_catmull_clark(
        &mut self,
        iterations: usize,
        update_points: bool,
    ) -> Result<(), Error> {
        if iterations == 0 {
            return Ok(());
        }
        check_for_deleted(&self.topol)?;
        let mut fpos = Vec::new();
        let mut epos = Vec::new();
        let mut vpos = if update_points {
            Some(Vec::new())
        } else {
            None
        };
        // Use vectors instead of properties because we don't want these to
        // change when we add new topology.
        CatmullClark::<DIM, A>::reserve(
            iterations,
            &mut self.topol,
            &mut vpos,
            &mut epos,
            &mut fpos,
        )?;
        // Temporary storage to use inside the loop.
        let mut hloop = Vec::new();
        let mut spliths = Vec::new();
        let mut subfaces = Vec::new();
        for _ in 0..iterations {
            CatmullClark::<DIM, A>::calc_face_edge_points(
                self,
                update_points,
                &mut fpos,
                &mut epos,
            )?;
            // Make them immutable from here.
            let fpos: &[A::Vector] = &fpos;
            let epos: &[A::Vector] = &epos;
            if let Some(vpos) = &mut vpos {
                CatmullClark::<DIM, A>::update_vertex_positions(self, fpos, epos, vpos)?;
            }
            let num_old_verts = self.num_vertices() as u32;
            for (ei, pos) in epos.iter().enumerate() {
                let ev = self.add_vertex(*pos)?;
                self.split_edge((ei as u32).into(), ev, true)?;
            }
            for f in self.faces() {
                CatmullClark::<DIM, A>::split_face(
                    self,
                    f,
                    num_old_verts,
                    fpos,
                    &mut hloop,
                    &mut spliths,
                    &mut subfaces,
                )?;
            }
        }
        Ok(())
    }
}

mod loop_scheme {
    use crate::{
        subdiv::check_for_deleted, topol::Topology, EditableTopology, Error, FloatScalarAdaptor,
        Handle, HasIterators, HasTopology, PolyMeshT,
    };
    use core::f64;
    use std::{
        marker::PhantomData,
        ops::{Add, Div, Mul},
    };

    const WEIGHTS: [(f64, f64); 256] = [
        (f64::MAX, f64::MAX),
        (0.765625, 0.234375),
        (0.390625, 0.3046875),
        (0.4375, 0.1875),
        (0.515625, 0.12109375),
        (0.5795339053710855, 0.08409321892578289),
        (0.625, 0.0625),
        (0.6568255586623777, 0.0490249201910889),
        (0.6794575214724776, 0.0400678098159403),
        (0.6959348386369, 0.033785017929233344),
        (0.7082224675195198, 0.02917775324804802),
        (0.7175917565621553, 0.025673476676167695),
        (0.7248797632095823, 0.022926686399201476),
        (0.7306500281453254, 0.020719228604205737),
        (0.735290719039789, 0.018907805782872215),
        (0.7390751047567019, 0.01739499301621987),
        (0.7421994992579459, 0.01611253129637838),
        (0.7448075716764623, 0.015011319313149278),
        (0.7470062552448259, 0.014055208041954117),
        (0.7488763737375063, 0.01321703296118388),
        (0.7504798778795585, 0.012476006106022078),
        (0.7518648627822762, 0.011815958915129706),
        (0.7530691054536927, 0.011224131570286695),
        (0.7541225977697901, 0.010690321836096082),
        (0.755049386297464, 0.010206275570939),
        (0.755868926462989, 0.009765242941480439),
        (0.7565970915690476, 0.0093616503242674),
        (0.7572469332438236, 0.008990854324302829),
        (0.7578292606560425, 0.00864895497656991),
        (0.7583530860603586, 0.008332652204815221),
        (0.7588259706889198, 0.00803913431036934),
        (0.7592542955979749, 0.007765990464581454),
        (0.7596434754665834, 0.007511141391669269),
        (0.7599981286435098, 0.0072727839804997045),
        (0.7603222133596177, 0.007049346665893595),
        (0.7606191375697349, 0.006839453212293289),
        (0.7608918480893487, 0.006641893108629203),
        (0.7611429033608542, 0.0064555972064634),
        (0.7613745331911552, 0.00627961754760118),
        (0.761588688055457, 0.006113110562680592),
        (0.761787079995812, 0.005955323000104702),
        (0.7619712167105319, 0.005805580080230927),
        (0.7621424300980284, 0.005663275473856465),
        (0.7623019002612208, 0.005527862784622771),
        (0.762450675778128, 0.005398848277769819),
        (0.7625896908871168, 0.005275784646952959),
        (0.7627197801114307, 0.0051582656497515065),
        (0.7628416907494706, 0.005045921473415519),
        (0.7629560935791229, 0.004938414717101607),
        (0.7630635920618283, 0.004835436896697382),
        (0.7631647302817344, 0.004736705394365313),
        (0.763259999814582, 0.004641960787949372),
        (0.7633498456879493, 0.004550964506000977),
        (0.7634346715675473, 0.004463496762876466),
        (0.7635148442822338, 0.004379354735514189),
        (0.763590697782301, 0.0042983509494127084),
        (0.7636625366106651, 0.004220311846238125),
        (0.7637306389542229, 0.004145076509575037),
        (0.7637952593323827, 0.004072495528752024),
        (0.763856630971209, 0.004002429983538829),
        (0.7639149679044828, 0.003934750534925288),
        (0.7639704668369653, 0.0038693366092300754),
        (0.7640233088001207, 0.0038060756645141814),
        (0.7640736606262798, 0.0037448625297415907),
        (0.764121676263638, 0.0036855988083806567),
        (0.7641674979514157, 0.003628192339208989),
        (0.7642112572719131, 0.003572556708001318),
        (0.7642530760939676, 0.003518610804567648),
        (0.7642930674204409, 0.0034662784202876347),
        (0.7643313361507167, 0.0034154878818736718),
        (0.7643679797678146, 0.003366171717602649),
        (0.7644030889585124, 0.0033182663526970086),
        (0.7644367481738337, 0.003271711830918976),
        (0.7644690361363657, 0.003226451559775811),
        (0.7645000263000896, 0.0031824320770258165),
        (0.7645297872677354, 0.0031396028364301955),
        (0.7645583831700856, 0.003097916010919926),
        (0.764585874011139, 0.00305732631154365),
        (0.7646123159825993, 0.0030177908207359074),
        (0.7646377617507658, 0.0029792688385979005),
        (0.7646622607185596, 0.002941721741018005),
        (0.7646858592651113, 0.002905112848578872),
        (0.7647086009650872, 0.002869407305303815),
        (0.7647305267896813, 0.0028345719663893814),
        (0.7647516752910065, 0.002800575294154684),
        (0.7647720827714336, 0.002767387261512546),
        (0.764791783439263, 0.0027349792623341503),
        (0.7648108095519792, 0.0027033240281381693),
        (0.7648291915482017, 0.002672395550588617),
        (0.7648469581693449, 0.0026421690093332032),
        (0.7648641365718911, 0.0026126207047567646),
        (0.7648807524310953, 0.0025837279952626887),
        (0.7648968300368615, 0.002555469238729767),
        (0.7649123923824582, 0.00252782373782303),
        (0.7649274612466781, 0.002500771688865126),
        (0.7649420572699882, 0.002474294134000124),
        (0.7649562000251697, 0.0024483729164044813),
        (0.7649699080828968, 0.002422990638320652),
        (0.7649831990726645, 0.0023981306217075054),
        (0.7649960897394407, 0.002373776871318781),
        (0.7650085959963784, 0.0023499140400362155),
        (0.7650207329739014, 0.0023265273962980062),
        (0.7650325150654438, 0.0023036027934760415),
        (0.7650439559701021, 0.002281126641066969),
        (0.7650550687324362, 0.002259085877572728),
        (0.7650658657796328, 0.002237467944955879),
        (0.7650763589562327, 0.002216260764563842),
        (0.7650865595565998, 0.002195452714424301),
        (0.7650964783552985, 0.0021750326078213097),
        (0.765106125635534, 0.0021549896730684946),
        (0.7651155112157926, 0.002135313534401886),
        (0.7651246444748127, 0.0021159941939206065),
        (0.7651335343750045, 0.0020970220145088885),
        (0.7651421894844288, 0.0020783877036776206),
        (0.7651506179974327, 0.0020600822982681346),
        (0.7651588277540367, 0.002042097149964899),
        (0.7651668262581595, 0.0020244239115675914),
        (0.7651746206947587, 0.0020070545239764216),
        (0.7651822179459596, 0.0019899812038478),
        (0.7651896246062413, 0.0019731964318803253),
        (0.7651968469967412, 0.001956692941693824),
        (0.7652038911787375, 0.0019404637092666319),
        (0.7652107629663616, 0.001924501942898675),
        (0.765217467938591, 0.0019088010736699912),
        (0.7652240114505691, 0.0018933547463663785),
        (0.7652303986442939, 0.0018781568108456486),
        (0.7652366344587167, 0.0018632013138197088),
        (0.7652427236392867, 0.001848482491029239),
        (0.7652486707469761, 0.0018339947597892496),
        (0.765254480166818, 0.0018197327118851326),
        (0.7652601561159877, 0.0018056911068000947),
        (0.7652657026514534, 0.0017918648652560805),
        (0.7652711236772229, 0.0017782490630513413),
        (0.7652764229512121, 0.0017648389251788562),
        (0.7652816040917552, 0.0017516298202107823),
        (0.7652866705837797, 0.0017386172549349644),
        (0.7652916257846647, 0.0017257968692304065),
        (0.7652964729298011, 0.0017131644311693356),
        (0.7653012151378717, 0.0017007158323342627),
        (0.7653058554158665, 0.0016884470833390896),
        (0.7653103966638477, 0.0016763543095439446),
        (0.7653148416794797, 0.0016644337469540445),
        (0.765319193162337, 0.0016526817382934015),
        (0.7653234537180025, 0.001641094729244738),
        (0.7653276258619655, 0.0016296692648474626),
        (0.765331712023335, 0.0016184019860459651),
        (0.7653357145483747, 0.0016072896263809953),
        (0.7653396357038686, 0.0015963290088172205),
        (0.7653434776803318, 0.0015855170427004608),
        (0.7653472425950671, 0.001574850720838476),
        (0.7653509324950829, 0.0015643271166994469),
        (0.7653545493598743, 0.0015539433817226868),
        (0.7653580951040773, 0.0015436967427363334),
        (0.7653615715800024, 0.001533584499477109),
        (0.7653649805800529, 0.0015236040222074492),
        (0.7653683238390346, 0.001513752749425583),
        (0.7653716030363633, 0.0015040281856643373),
        (0.7653748197981722, 0.0014944278993746988),
        (0.7653779756993283, 0.0014849495208903271),
        (0.765381072265359, 0.0014755907404694403),
        (0.7653841109742958, 0.0014663493064106512),
        (0.7653870932584372, 0.00145722302323952),
        (0.7653900205060367, 0.0014482097499627367),
        (0.7653928940629173, 0.0014393073983870104),
        (0.7653957152340187, 0.001430513931499886),
        (0.7653984852848787, 0.0014218273619098257),
        (0.7654012054430512, 0.0014132457503430648),
        (0.7654038768994671, 0.0014047672041948079),
        (0.7654065008097366, 0.0013963898761325203),
        (0.7654090782953994, 0.001388111962749116),
        (0.7654116104451227, 0.0013799317032639842),
        (0.7654140983158503, 0.0013718473782698814),
        (0.7654165429339058, 0.0013638573085238033),
        (0.7654189452960509, 0.0013559598537800528),
        (0.7654213063705004, 0.0013481534116637909),
        (0.7654236270978987, 0.0013404364165834357),
        (0.7654259083922563, 0.0013328073386803614),
        (0.7654281511418488, 0.0013252646828144136),
        (0.7654303562100823, 0.0013178069875838074),
        (0.7654325244363238, 0.0013104328243780795),
        (0.7654346566367001, 0.0013031407964627776),
        (0.7654367536048652, 0.0012959295380946673),
        (0.7654388161127397, 0.0012887977136662651),
        (0.7654408449112199, 0.00128174401687858),
        (0.765442840730862, 0.0012747671699409673),
        (0.765444804282539, 0.0012678659227970863),
        (0.7654467362580744, 0.0012610390523759446),
        (0.7654486373308506, 0.0012542853618671094),
        (0.7654505081563963, 0.001247603680019169),
        (0.7654523493729513, 0.0012409928604605756),
        (0.7654541616020104, 0.0012344517810420504),
        (0.7654559454488484, 0.0012279793431997465),
        (0.7654577015030251, 0.0012215744713384115),
        (0.7654594303388719, 0.0012152361122338244),
        (0.7654611325159626, 0.001208963234453801),
        (0.7654628085795647, 0.001202754827797104),
        (0.7654644590610769, 0.0011966099027496077),
        (0.7654660844784502, 0.0011905274899571053),
        (0.7654676853365936, 0.0011845066397141737),
        (0.7654692621277668, 0.0011785464214685085),
        (0.7654708153319583, 0.001172645923340209),
        (0.7654723454172504, 0.0011668042516554706),
        (0.7654738528401726, 0.0011610205304941955),
        (0.7654753380460413, 0.0011552939012510284),
        (0.7654768014692888, 0.0011496235222093689),
        (0.7654782435337811, 0.001144008568127897),
        (0.7654796646531251, 0.0011384482298391985),
        (0.7654810652309644, 0.0011329417138600756),
        (0.7654824456612661, 0.0011274882420131435),
        (0.7654838063285988, 0.001122087051059336),
        (0.7654851476083986, 0.0011167373923409594),
        (0.7654864698672297, 0.0011114385314349307),
        (0.7654877734630341, 0.0011061897478158772),
        (0.7654890587453733, 0.001100990334528764),
        (0.765490326055664, 0.0010958395978707285),
        (0.7654915757274035, 0.0010907368570818442),
        (0.7654928080863899, 0.0010856814440444917),
        (0.7654940234509333, 0.0010806727029910907),
        (0.7654952221320632, 0.0010757099902198935),
        (0.7654964044337258, 0.0010707926738186037),
        (0.7654975706529772, 0.001065920133395558),
        (0.7654987210801703, 0.0010610917598182338),
        (0.7654998559991347, 0.0010563069549588526),
        (0.765500975687353, 0.0010515651314468477),
        (0.7655020804161292, 0.0010468657124279945),
        (0.7655031704507533, 0.0010422081313299853),
        (0.7655042460506614, 0.001037591831634242),
        (0.7655053074695883, 0.0010330162666537961),
        (0.7655063549557195, 0.0010284808993170197),
        (0.7655073887518341, 0.0010239852019570564),
        (0.7655084090954467, 0.0010195286561067538),
        (0.7655094162189431, 0.0010151107522989474),
        (0.7655104103497133, 0.0010107309898719251),
        (0.7655113917102796, 0.001006388876779916),
        (0.7655123605184204, 0.0010020839294084597),
        (0.7655133169872923, 0.0009978156723945005),
        (0.7655142613255461, 0.0009935836384510758),
        (0.7655151937374416, 0.000989387368196449),
        (0.7655161144229572, 0.000985226409987575),
        (0.7655170235778974, 0.0009811003197577515),
        (0.7655179213939973, 0.0009770086608583443),
        (0.7655188080590236, 0.0009729510039044664),
        (0.7655196837568718, 0.0009689269266244965),
        (0.7655205486676627, 0.0009649360137133222),
        (0.7655214029678348, 0.000960977856689202),
        (0.7655222468302332, 0.0009570520537541501),
        (0.765523080424199, 0.0009531582096577273),
        (0.7655239039156523, 0.0009492959355641606),
        (0.7655247174671761, 0.0009454648489226771),
        (0.7655255212380956, 0.0009416645733409814),
        (0.7655263153845568, 0.0009378947384617728),
        (0.7655271000596022, 0.0009341549798422224),
        (0.7655278754132444, 0.0009304449388363319),
        (0.765528641592538, 0.0009267642624800873),
        (0.7655293987416487, 0.0009231126033793358),
        (0.7655301470019225, 0.0009194896196003039),
    ];

    pub fn compute_weight(valence: usize) -> (f64, f64) {
        let alpha = (40.0
            - f64::powf(
                3.0 + 2.0 * f64::cos(f64::consts::TAU / (valence as f64)),
                2.0,
            ))
            / 64.0;
        (1.0 - alpha, alpha / (valence as f64))
    }

    /// This struct doesn't conatin any data. This purpose of this struct is to
    /// provide a convenient scope inside its `impl` where I can impose all the
    /// trait bounds once and write several free functions within that scope, with
    /// those trait bounds applied.
    struct LoopScheme<const DIM: usize, A>(PhantomData<A>);

    impl<const DIM: usize, A> LoopScheme<DIM, A>
    where
        A: FloatScalarAdaptor<DIM>,
        A::Scalar: Add<Output = A::Scalar>,
        A::Vector: Add<Output = A::Vector>
            + Div<A::Scalar, Output = A::Vector>
            + Mul<A::Scalar, Output = A::Vector>,
    {
        fn reserve(
            iterations: usize,
            mesh: &mut Topology,
            vertex_points: &mut Vec<A::Vector>,
            edge_points: &mut Vec<A::Vector>,
        ) -> Result<(), Error> {
            let (mut nv, mut ne, mut nf) =
                (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces());
            for i in 0..iterations {
                if i == iterations - 1 {
                    vertex_points.reserve(nv);
                    edge_points.reserve(ne);
                }
                let v = nv + ne;
                let e = ne * 2 + nf * 3;
                let f = nf * 4;
                (nv, ne, nf) = (v, e, f);
            }
            mesh.reserve(nv, ne, nf)
        }

        fn compute_vertex_points(mesh: &Topology, points: &[A::Vector], dst: &mut Vec<A::Vector>) {
            dst.clear();
            dst.extend(mesh.vertices().map(|v| {
                if v.is_boundary(mesh) {
                    if let Some(h) = v.halfedge(mesh) {
                        debug_assert!(h.is_boundary(mesh));
                        (points[v.index() as usize] * A::scalarf64(6.0)
                            + points[h.head(mesh).index() as usize]
                            + points[h.prev(mesh).tail(mesh).index() as usize])
                            / A::scalarf64(8.0)
                    } else {
                        // Isolated vertex doesn't move.
                        points[v.index() as usize]
                    }
                } else {
                    let (valence, sum) = mesh
                        .vv_cw_iter(v)
                        .fold((0usize, A::zero_vector()), |(valence, sum), nv| {
                            (valence + 1, sum + points[nv.index() as usize])
                        });
                    let (a, b) = if valence < WEIGHTS.len() {
                        WEIGHTS[valence]
                    } else {
                        compute_weight(valence)
                    };
                    sum * A::scalarf64(b) + points[v.index() as usize] * A::scalarf64(a)
                }
            }));
        }
    }

    impl<const DIM: usize, A> PolyMeshT<DIM, A>
    where
        A: FloatScalarAdaptor<DIM>,
        A::Scalar: Add<Output = A::Scalar>,
        A::Vector: Add<Output = A::Vector>
            + Div<A::Scalar, Output = A::Vector>
            + Mul<A::Scalar, Output = A::Vector>,
    {
        pub fn subdivide_loop(
            &mut self,
            iterations: usize,
            update_points: bool,
        ) -> Result<(), Error> {
            if iterations == 0 {
                return Ok(());
            }
            check_for_deleted(&self.topol)?;
            self.triangulate()?;
            let mut vpos = Vec::new();
            let mut epos = Vec::new();
            LoopScheme::<DIM, A>::reserve(iterations, &mut self.topol, &mut vpos, &mut epos)?;
            let mut hhs = Vec::new();
            for _ in 0..iterations {
                {
                    let mut points = self.points();
                    let mut points = points.try_borrow_mut()?;
                    // Compute vertex points.
                    if update_points {
                        LoopScheme::<DIM, A>::compute_vertex_points(
                            &self.topol,
                            &points,
                            &mut vpos,
                        );
                        points.copy_from_slice(&vpos);
                    }
                }
                {
                    let points = self.points();
                    let points = points.try_borrow()?;
                    // Compute edge points.
                    epos.clear();
                    epos.extend(self.edges().map(|e| {
                        if e.is_boundary(self) || !update_points {
                            let (v0, v1) = e.vertices(self);
                            (points[v0.index() as usize] + points[v1.index() as usize])
                                * A::scalarf64(0.5)
                        } else {
                            let (h, oh) = e.halfedges();
                            ((points[h.head(self).index() as usize]
                                + points[oh.head(self).index() as usize])
                                * A::scalarf64(3.0)
                                + points[h.next(self).head(self).index() as usize]
                                + points[oh.next(self).head(self).index() as usize])
                                / A::scalarf64(8.0)
                        }
                    }));
                }
                // Make them immutable from here.
                let epos: &[A::Vector] = &epos;
                let num_old_verts = self.num_vertices() as u32;
                // Split edges.
                for (e, epos) in self.edges().map(|e| (e, epos[e.index() as usize])) {
                    let ev = self.add_vertex(epos)?;
                    self.split_edge(e, ev, true)?;
                }
                for f in self.faces() {
                    let hstart = self
                        .fh_ccw_iter(f)
                        .find(|h| h.head(self).index() < num_old_verts)
                        .ok_or(Error::CannotSplitFace(f))?;
                    hhs.clear();
                    hhs.extend(self.loop_ccw_iter(hstart));
                    debug_assert!(hhs.len() % 2 == 0);
                    for hpair in hhs.chunks_exact(2) {
                        self.insert_edge(hpair[1], hpair[0])?;
                        self.check_topology().expect("Failed...");
                    }
                }
            }
            Ok(())
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{alum_glam::PolyMeshF32, obj::test::bunny_mesh, HasTopology};

    #[test]
    fn t_box_catmull_clark() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
        mesh.subdivide_catmull_clark(1, true)
            .expect("Subdivision failed");
        assert_eq!(26, mesh.num_vertices());
        assert_eq!(48, mesh.num_edges());
        assert_eq!(24, mesh.num_faces());
        mesh.check_topology().expect("Topological errors found");
    }

    #[test]
    fn t_bunny_subdiv_catmull_clark() {
        let mut mesh = bunny_mesh();
        mesh.subdivide_catmull_clark(3, true)
            .expect("Cannot subivide");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(mesh.try_calc_area().expect("Cannot compute area"), 5.566642);
        assert_eq!(238464, mesh.num_faces());
        assert_eq!(477096, mesh.num_edges());
        assert_eq!(238630, mesh.num_vertices());
    }

    #[test]
    fn t_triangle_subdiv_loop() {
        let mut mesh = PolyMeshF32::new();
        mesh.add_vertex(glam::vec3(0.0, 0.0, 0.0))
            .expect("Cannot add vertex");
        mesh.add_vertex(glam::vec3(1.0, 0.0, 0.0))
            .expect("Cannot add vertex");
        mesh.add_vertex(glam::vec3(1.0, 1.0, 0.0))
            .expect("Cannot add vertex");
        mesh.add_tri_face(0.into(), 1.into(), 2.into())
            .expect("Cannot add face");
        mesh.subdivide_loop(1, true).expect("Cannot subidivde");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(6, mesh.num_vertices());
        assert_eq!(9, mesh.num_edges());
        assert_eq!(4, mesh.num_faces());
    }

    #[test]
    fn t_box_subdiv_loop() {
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot create a box");
        mesh.subdivide_loop(1, true).expect("Cannot subdivide");
        mesh.check_topology().expect("Topological errors found");
        // assert_eq!(mesh.try_calc_area().expect("Cannot compute area"), 5.566642);
        // DEBUG
        for vi in 0..8 {
            eprintln!("{}: {}", vi, mesh.point(vi.into()).unwrap());
        }
        todo!()
    }
}
