use super::DecimaterModule;
use crate::{
    CrossProductAdaptor, Error, FloatScalarAdaptor, Handle, HasIterators, PolyMeshT,
    VectorNormalizeAdaptor, HH, VH,
};
use std::{
    marker::PhantomData,
    ops::{Add, Div, Mul, Sub},
};

pub struct QuadricModule<A>(PhantomData<A>)
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Add<Output = A::Vector> + Sub<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>;

impl<A> DecimaterModule<PolyMeshT<3, A>> for QuadricModule<A>
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Add<Output = A::Vector> + Sub<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    fn collapse_cost(&self, _mesh: &PolyMeshT<3, A>, _h: HH) -> Option<f64> {
        todo!()
    }

    fn before_collapse(&self, _mesh: &PolyMeshT<3, A>, _h: HH) -> Result<(), Error> {
        todo!()
    }

    fn after_collapse(&self, mesh: &mut PolyMeshT<3, A>, v: VH) -> Result<(), Error> {
        if let Some(mut fnormals) = mesh.face_normals() {
            let mut fnormals = fnormals.try_borrow_mut()?;
            let points = mesh.points();
            let points = points.try_borrow()?;
            for f in mesh.vf_ccw_iter(v) {
                fnormals[f.index() as usize] = mesh.calc_face_normal(f, &points);
            }
            if let Some(mut vnormals) = mesh.vertex_normals() {
                let mut vnormals = vnormals.try_borrow_mut()?;
                vnormals[v.index() as usize] = mesh.calc_vertex_normal_fast(v, &fnormals);
            }
        }
        todo!()
    }
}
