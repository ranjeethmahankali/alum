use core::f64;
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul},
};

use crate::{
    topol::Topology, EditableTopology, Error, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, PolyMeshT,
};

use super::check_for_deleted;

const WEIGHTS: [(f64, f64); 64] = [
    (f64::MAX, f64::MAX), // Should never happen.
    (0.765625, 0.234375),
    (0.390625, 0.3046875),
    (0.4375, 0.1875),
    (0.515625, 0.12109375),
    (0.5795339053710855, 0.08409321892578289),
    (0.625, 0.0625),
    (0.6568255586623777, 0.0490249201910889),
    (0.6794575214724776, 0.0400678098159403),
    (0.6959348386369, 0.03378501792923334),
    (0.7082224675195198, 0.02917775324804802),
    (0.7175917565621553, 0.025673476676167695),
    (0.7248797632095823, 0.02292668639920148),
    (0.7306500281453254, 0.020719228604205737),
    (0.735290719039789, 0.01890780578287221),
    (0.7390751047567019, 0.017394993016219874),
    (0.7421994992579459, 0.01611253129637838),
    (0.7448075716764623, 0.015011319313149278),
    (0.7470062552448259, 0.014055208041954117),
    (0.7488763737375063, 0.01321703296118388),
    (0.7504798778795585, 0.012476006106022074),
    (0.7518648627822762, 0.011815958915129706),
    (0.7530691054536927, 0.011224131570286698),
    (0.7541225977697901, 0.010690321836096084),
    (0.755049386297464, 0.010206275570939),
    (0.755868926462989, 0.00976524294148044),
    (0.7565970915690476, 0.009361650324267399),
    (0.7572469332438236, 0.008990854324302829),
    (0.7578292606560425, 0.00864895497656991),
    (0.7583530860603586, 0.008332652204815221),
    (0.7588259706889198, 0.008039134310369341),
    (0.7592542955979749, 0.007765990464581456),
    (0.7596434754665834, 0.007511141391669269),
    (0.7599981286435098, 0.007272783980499703),
    (0.7603222133596177, 0.007049346665893597),
    (0.7606191375697349, 0.006839453212293289),
    (0.7608918480893487, 0.006641893108629202),
    (0.7611429033608542, 0.0064555972064634),
    (0.7613745331911552, 0.006279617547601178),
    (0.761588688055457, 0.00611311056268059),
    (0.761787079995812, 0.005955323000104701),
    (0.7619712167105319, 0.005805580080230929),
    (0.7621424300980284, 0.005663275473856467),
    (0.7623019002612208, 0.005527862784622772),
    (0.762450675778128, 0.005398848277769819),
    (0.7625896908871168, 0.0052757846469529595),
    (0.7627197801114307, 0.0051582656497515065),
    (0.7628416907494706, 0.005045921473415519),
    (0.7629560935791229, 0.004938414717101607),
    (0.7630635920618283, 0.004835436896697382),
    (0.7631647302817344, 0.004736705394365313),
    (0.763259999814582, 0.004641960787949372),
    (0.7633498456879493, 0.004550964506000976),
    (0.7634346715675473, 0.004463496762876466),
    (0.7635148442822338, 0.004379354735514189),
    (0.763590697782301, 0.0042983509494127084),
    (0.7636625366106651, 0.004220311846238124),
    (0.7637306389542229, 0.004145076509575037),
    (0.7637952593323827, 0.004072495528752023),
    (0.763856630971209, 0.004002429983538829),
    (0.7639149679044828, 0.003934750534925287),
    (0.7639704668369653, 0.0038693366092300762),
    (0.7640233088001207, 0.0038060756645141823),
    (0.7640736606262798, 0.0037448625297415907),
];

fn compute_weights(valence: usize) -> (f64, f64) {
    const THREE_OVER_EIGHT: f64 = 3.0 / 8.0;
    let n = valence as f64;
    let alpha =
        THREE_OVER_EIGHT + f64::powi(THREE_OVER_EIGHT + 0.25 * f64::cos(f64::consts::TAU / n), 2);
    (alpha, (1.0 - alpha) / n)
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
                let (a, b) = match WEIGHTS.get(valence) {
                    Some((a, b)) => (*a, *b),
                    None => compute_weights(valence),
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
    /// Subdivide the mesh according to teh [Loop subdivision
    /// scheme](https://en.wikipedia.org/wiki/Loop_subdivision_surface).
    ///
    /// Subdivisions are carried out for the given number of
    /// `iterations`. `update_points` determines whether the vertices of the
    /// mesh are moved. If this is `false`, only the topology of th emesh is
    /// subdivided, leaving the shape of the mesh unmodified. If this is `true`
    /// then the shape is smoothed according to the Loop-Subdivision
    /// scheme. Loop subdivision is meant for triangle meshes. If the mesh is
    /// not a triangle mesh, it will be triangulated before any subdivision is
    /// performed.
    ///
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, HasTopology};
    ///
    /// let mut mesh = PolyMeshF32::unit_box().expect("Cannot create box");
    /// assert_eq!((8, 12, 6), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.subdivide_loop(1, true).expect("Subdivision failed");
    /// // The mesh now has more faces.
    /// assert_eq!((26, 72, 48), (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces()));
    /// mesh.check_topology().expect("Topological errors found");
    /// ```
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
                points.copy_from_slice(vpos);
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
