use crate::{Adaptor, Error, FloatScalarAdaptor, Handle, PolyMeshT};
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul},
};

/// This struct doesn't conatin any data. This purpose of this struct is to
/// provide a convenient scope inside its `impl` where I can impose all the
/// trait bounds once and write several free functions within that scope, with
/// those trait bounds applied.
struct CatmullClark<const DIM: usize, A>(PhantomData<A>);

impl<const DIM: usize, A> CatmullClark<DIM, A>
where
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    fn calc_face_edge_points(
        mesh: &mut PolyMeshT<DIM, A>,
        update_points: bool,
        fpos: &mut Vec<A::Vector>,
        epos: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        // Compute face points.
        fpos.clear();
        fpos.reserve(mesh.num_faces());
        fpos.extend(mesh.faces().map(|f| mesh.calc_face_centroid(f, &points)));
        // Compute edge points.
        epos.clear();
        epos.reserve(mesh.num_edges());
        epos.extend(mesh.edges().map(|e| {
            let (h, oh) = mesh.halfedge_pair(e);
            match (mesh.halfedge_face(h), mesh.halfedge_face(oh)) {
                (Some(fa), Some(fb)) if update_points => {
                    (points[mesh.head_vertex(h).index() as usize]
                        + points[mesh.head_vertex(oh).index() as usize]
                        + fpos[fa.index() as usize]
                        + fpos[fb.index() as usize])
                        * A::scalarf64(0.25)
                }
                _ => {
                    (points[mesh.head_vertex(h).index() as usize]
                        + points[mesh.head_vertex(oh).index() as usize])
                        * A::scalarf64(0.5)
                }
            }
        }));
        Ok(())
    }

    fn update_vertex_positions(
        mesh: &mut PolyMeshT<DIM, A>,
        fpos: &[A::Vector],
        epos: &[A::Vector],
        vpos: &mut Vec<A::Vector>,
    ) -> Result<(), Error> {
        // Compute vertex positions.
        let points = mesh.points();
        let points = points.try_borrow()?;
        vpos.clear();
        vpos.reserve(mesh.num_vertices());
        vpos.extend(mesh.vertices().map(|v| {
            if mesh.is_boundary_vertex(v) {
                let (count, sum) = mesh
                    .ve_ccw_iter(v)
                    .filter(|e| mesh.is_boundary_edge(*e))
                    .fold((1usize, points[v.index() as usize]), |(count, total), e| {
                        (count + 1, total + epos[e.index() as usize])
                    });
                sum / A::scalarf64(count as f64)
            } else {
                let valence = mesh.vertex_valence(v) as f64;
                (((mesh.vf_ccw_iter(v).fold(A::zero_vector(), |total, f| {
                    total + fpos[f.index() as usize]
                }) + mesh.vv_ccw_iter(v).fold(A::zero_vector(), |total, v| {
                    total + points[v.index() as usize]
                })) / A::scalarf64(valence))
                    + (points[v.index() as usize] * A::scalarf64(valence - 2.0)))
                    / A::scalarf64(valence)
            }
        }));
        // Update the vertex positions.
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut()?;
        points.copy_from_slice(vpos);
        Ok(())
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    pub fn subidivide_catmull_clark(
        &mut self,
        iters: usize,
        update_points: bool,
    ) -> Result<(), Error> {
        if iters == 0 {
            return Ok(());
        }
        {
            // Ensure no deleted topology.
            let fstatus = self.topol.fstatus.try_borrow()?;
            if let Some(f) = self.faces().find(|f| fstatus[f.index() as usize].deleted()) {
                return Err(Error::DeletedFace(f));
            }
            let estatus = self.topol.estatus.try_borrow()?;
            if let Some(e) = self.edges().find(|e| estatus[e.index() as usize].deleted()) {
                return Err(Error::DeletedEdge(e));
            }
            let vstatus = self.topol.hstatus.try_borrow()?;
            if let Some(v) = self
                .vertices()
                .find(|v| vstatus[v.index() as usize].deleted())
            {
                return Err(Error::DeletedVertex(v));
            }
        }
        // Use vectors instead of properties because we don't want these to
        // change when we add new topology.
        let mut fpos = Vec::with_capacity(self.num_faces());
        let mut epos = Vec::with_capacity(self.num_edges());
        let mut vpos = Vec::new();
        for _ in 0..iters {
            CatmullClark::<DIM, A>::calc_face_edge_points(
                self,
                update_points,
                &mut fpos,
                &mut epos,
            )?;
            if update_points {
                CatmullClark::<DIM, A>::update_vertex_positions(self, &fpos, &epos, &mut vpos)?;
            }
            todo!("Modify the topolgy");
        }
        todo!()
    }
}
