use crate::mesh::{
    self, Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, VectorAngleAdaptor,
    VectorLengthAdaptor, VectorNormalizeAdaptor,
};

pub struct BuiltInAdaptorF32 {}

impl Adaptor<3> for BuiltInAdaptorF32 {
    type Vector = glam::Vec3;
    type Scalar = f32;

    fn vector(coords: [Self::Scalar; 3]) -> Self::Vector {
        glam::vec3(coords[0], coords[1], coords[2])
    }

    fn zero_vector() -> Self::Vector {
        glam::Vec3::splat(0.)
    }

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
    }
}

impl VectorLengthAdaptor<3> for BuiltInAdaptorF32 {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        glam::Vec3::length(v)
    }
}

impl VectorNormalizeAdaptor<3> for BuiltInAdaptorF32 {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for BuiltInAdaptorF32 {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(b)
    }
}

impl VectorAngleAdaptor for BuiltInAdaptorF32 {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl CrossProductAdaptor for BuiltInAdaptorF32 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

impl FloatScalarAdaptor<3> for BuiltInAdaptorF32 {
    fn scalarf32(val: f32) -> Self::Scalar {
        val
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val as f32
    }
}

pub struct BuiltInAdaptorF64 {}

impl Adaptor<3> for BuiltInAdaptorF64 {
    type Vector = glam::DVec3;
    type Scalar = f64;

    fn vector(coords: [Self::Scalar; 3]) -> Self::Vector {
        glam::dvec3(coords[0], coords[1], coords[2])
    }

    fn zero_vector() -> Self::Vector {
        glam::DVec3::splat(0.)
    }

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
    }
}

impl VectorLengthAdaptor<3> for BuiltInAdaptorF64 {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        glam::DVec3::length(v)
    }
}

impl VectorNormalizeAdaptor<3> for BuiltInAdaptorF64 {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for BuiltInAdaptorF64 {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(b)
    }
}

impl VectorAngleAdaptor for BuiltInAdaptorF64 {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl CrossProductAdaptor for BuiltInAdaptorF64 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

impl FloatScalarAdaptor<3> for BuiltInAdaptorF64 {
    fn scalarf32(val: f32) -> Self::Scalar {
        val as f64
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val
    }
}

pub type PolyMeshF32 = mesh::PolyMeshT<3, BuiltInAdaptorF32>;
pub type PolyMeshF64 = mesh::PolyMeshT<3, BuiltInAdaptorF64>;
