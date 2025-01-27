use super::Decimater;
use crate::{Adaptor, Error, HasTopology, PolyMeshT, FH, HH, VH};

pub struct HausdorffDecimater<A>
where
    A: Adaptor<3>,
{
    ptbuf: Vec<A::Vector>,
    fpoints: Vec<Vec<A::Vector>>,
    max_tol: A::Scalar,
}

impl<A> HausdorffDecimater<A>
where
    A: Adaptor<3>,
{
    pub fn new(mesh: &PolyMeshT<3, A>, max_tol: A::Scalar) -> Self {
        HausdorffDecimater {
            ptbuf: Vec::with_capacity(mesh.num_vertices()),
            fpoints: vec![Vec::new(); mesh.num_faces()],
            max_tol,
        }
    }
}

impl<A> Decimater<PolyMeshT<3, A>> for HausdorffDecimater<A>
where
    A: Adaptor<3>,
    A::Scalar: PartialOrd,
{
    type Cost = A::Scalar;

    fn collapse_cost(&self, mesh: &PolyMeshT<3, A>, h: HH) -> Option<Self::Cost> {
        todo!()
    }

    fn before_collapse(&mut self, mesh: &PolyMeshT<3, A>, h: HH) -> Result<(), Error> {
        Ok(()) // Nothing to do.
    }

    fn after_collapse(&mut self, mesh: &PolyMeshT<3, A>, v: VH) -> Result<(), Error> {
        todo!()
    }
}

fn triangle_dist_sq<A>(p: &A::Vector, a: &A::Vector, b: &A::Vector, c: &A::Vector) -> A::Scalar
where
    A: Adaptor<3>,
{
    todo!("Implement point to triangle squared distance");
}
