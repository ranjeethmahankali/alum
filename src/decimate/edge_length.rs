use crate::{DotProductAdaptor, Error, PolyMeshT, HH, VH};
use std::ops::{Mul, Sub};

use super::Decimater;

pub struct EdgeLengthDecimater<const DIM: usize, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd,
{
    max_length: A::Scalar,
}

impl<const DIM: usize, A> EdgeLengthDecimater<DIM, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd + Mul<Output = A::Scalar>,
{
    pub fn new(max_length: A::Scalar) -> Self {
        EdgeLengthDecimater {
            max_length: max_length * max_length,
        }
    }
}

impl<const DIM: usize, A> Decimater<PolyMeshT<DIM, A>> for EdgeLengthDecimater<DIM, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: PartialOrd + Mul<Output = A::Scalar>,
{
    type Cost = A::Scalar;

    fn collapse_cost(&self, mesh: &PolyMeshT<DIM, A>, h: HH) -> Option<Self::Cost> {
        let points = mesh.points();
        let points = points.try_borrow().ok()?;
        let elen = mesh.calc_edge_length_sqr(h.edge(), &points);
        if elen < self.max_length {
            Some(elen)
        } else {
            None
        }
    }

    fn before_collapse(&mut self, _mesh: &PolyMeshT<DIM, A>, _h: HH) -> Result<(), Error> {
        Ok(()) // Do nothing.
    }

    fn after_collapse(&mut self, _mesh: &PolyMeshT<DIM, A>, _v: VH) -> Result<(), Error> {
        Ok(()) // Do nothing.
    }
}
