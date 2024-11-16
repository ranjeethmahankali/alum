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

const LOOP_WEIGHTS: [(f64, f64); 64] = [
    (f64::MAX, f64::MAX), // Should never happen.
    (0.234375, 0.765625),
    (0.609375, 0.1953125),
    (0.5625, 0.14583333333333334),
    (0.484375, 0.12890625),
    (0.4204660946289145, 0.1159067810742171),
    (0.375, 0.10416666666666667),
    (0.3431744413376223, 0.09383222266605395),
    (0.3205424785275224, 0.0849321901840597),
    (0.3040651613631, 0.07732609318187778),
    (0.2917775324804802, 0.07082224675195198),
    (0.28240824343784465, 0.06523561423292322),
    (0.27512023679041775, 0.06040664693413186),
    (0.2693499718546746, 0.05620384831887119),
    (0.26470928096021096, 0.05252076564569922),
    (0.2609248952432981, 0.049271673650446796),
    (0.2578005007420541, 0.04638746870362162),
    (0.2551924283235377, 0.043812210098615426),
    (0.2529937447551741, 0.04150034751360144),
    (0.25112362626249374, 0.03941454598618454),
    (0.2495201221204415, 0.037523993893977925),
    (0.24813513721772384, 0.035803088703917914),
    (0.24693089454630734, 0.03423041388425876),
    (0.24587740223020993, 0.03278793903346913),
    (0.244950613702536, 0.031460391095727666),
    (0.24413107353701102, 0.030234757058519558),
    (0.24340290843095236, 0.029099888137271063),
    (0.24275306675617636, 0.028046182712734208),
    (0.2421707393439575, 0.027065330737715804),
    (0.24164691393964144, 0.02615010641587443),
    (0.24117402931108023, 0.02529419902296399),
    (0.24074570440202514, 0.024492074051547578),
    (0.2403565245334166, 0.02373885860833073),
    (0.2400018713564902, 0.0230302463225306),
    (0.2396777866403823, 0.022362418039988757),
    (0.23938086243026513, 0.02173197535913528),
    (0.23910815191065127, 0.021135884669148577),
    (0.23885709663914578, 0.020571429820563626),
    (0.23862546680884478, 0.02003617192608303),
    (0.23841131194454301, 0.019527915078345053),
    (0.23821292000418803, 0.0190446769998953),
    (0.23802878328946808, 0.018584663822208095),
    (0.2378575699019716, 0.01814624833566734),
    (0.2376980997387792, 0.0177279511688656),
    (0.23754932422187203, 0.01732842444950291),
    (0.2374103091128832, 0.01694643757526926),
    (0.2372802198885693, 0.016580864785031103),
    (0.2371583092505294, 0.016230674271265334),
    (0.23704390642087714, 0.015894918616231726),
    (0.23693640793817172, 0.015572726368608741),
    (0.23683526971826563, 0.015263294605634687),
    (0.23674000018541796, 0.01496588234930553),
    (0.23665015431205072, 0.014679804724768255),
    (0.2365653284324527, 0.01440442776542542),
    (0.23648515571776618, 0.01413916378300433),
    (0.23640930221769896, 0.013883467232405473),
    (0.23633746338933492, 0.013636831010904734),
    (0.23626936104577712, 0.01339878313954777),
    (0.23620474066761732, 0.013168883781592805),
    (0.23614336902879096, 0.012946722558834051),
    (0.23608503209551723, 0.01273191613174138),
    (0.23602953316303465, 0.012524106013720743),
    (0.2359766911998793, 0.012322956593550333),
    (0.2359263393737202, 0.012128153343274283),
];

fn compute_loop_weights(valence: usize) -> (f64, f64) {
    const THREE_OVER_EIGHT: f64 = 3.0 / 8.0;
    let n = valence as f64;
    let alpha =
        THREE_OVER_EIGHT + f64::powi(THREE_OVER_EIGHT + 0.25 * f64::cos(f64::consts::TAU / n), 2);
    (1.0 - alpha, alpha / n)
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
        let (mut nv, mut ne, mut nf) = (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces());
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
                    .vv_ccw_iter(v)
                    .fold((0usize, A::zero_vector()), |(valence, sum), nv| {
                        (valence + 1, sum + points[nv.index() as usize])
                    });
                let (a, b) = match LOOP_WEIGHTS.get(valence) {
                    Some((a, b)) => (*a, *b),
                    None => compute_loop_weights(valence),
                };
                sum * A::scalarf64(b) + points[v.index() as usize] * A::scalarf64(a)
            }
        }));
    }

    fn compute_edge_points(
        mesh: &Topology,
        update_points: bool,
        points: &[A::Vector],
        dst: &mut Vec<A::Vector>,
    ) {
        dst.clear();
        dst.extend(mesh.edges().map(|e| {
            let (v0, v1) = e.vertices(mesh);
            let vsum = points[v0.index() as usize] + points[v1.index() as usize];
            if e.is_boundary(mesh) || !update_points {
                vsum * A::scalarf64(0.5)
            } else {
                let (h, oh) = e.halfedges();
                let (v2, v3) = (h.next(mesh).head(mesh), oh.next(mesh).head(mesh));
                (vsum * A::scalarf64(3.0)
                    + (points[v2.index() as usize] + points[v3.index() as usize]))
                    / A::scalarf64(8.0)
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
    pub fn subdivide_loop(&mut self, iterations: usize, update_points: bool) -> Result<(), Error> {
        if iterations == 0 {
            return Ok(());
        }
        check_for_deleted(&self.topol)?;
        self.triangulate()?;
        let mut vpos = Vec::new();
        let mut epos = Vec::new();
        LoopScheme::<DIM, A>::reserve(iterations, &mut self.topol, &mut vpos, &mut epos)?;
        let mut hhs = Vec::new();
        let mut points = self.points();
        for _ in 0..iterations {
            {
                let points = points.try_borrow()?;
                // Compute vertex points.
                if update_points {
                    LoopScheme::<DIM, A>::compute_vertex_points(&self.topol, &points, &mut vpos);
                }
                // Compute edge points.
                LoopScheme::<DIM, A>::compute_edge_points(
                    &self.topol,
                    update_points,
                    &points,
                    &mut epos,
                );
            }
            // Make them immutable from here.
            let epos: &[A::Vector] = &epos;
            let vpos: &[A::Vector] = &vpos;
            if update_points {
                let mut points = points.try_borrow_mut()?;
                points.copy_from_slice(&vpos);
            }
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
                }
            }
        }
        Ok(())
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
        let mut mesh = PolyMeshF32::unit_box().expect("Cannot subdivide");
        mesh.subdivide_loop(2, true).expect("Cannot subdivide");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(12 * 16, mesh.num_faces());
        assert_eq!(288, mesh.num_edges());
        assert_eq!(192, mesh.num_faces());
        assert_eq!(
            mesh.try_calc_area().expect("Cannot compute area"),
            2.6263301
        );
    }

    #[test]
    fn t_bunny_subdiv_loop() {
        let mut mesh = bunny_mesh();
        mesh.subdivide_loop(3, true)
            .expect("Cannot subdivide the mesh");
        assert_eq!(159142, mesh.num_vertices());
        assert_eq!(477096, mesh.num_edges());
        assert_eq!(317952, mesh.num_faces());
    }
}
