/*!
This is an optional module that is enabled by the `use_glam` feature. It
provides mesh types that can be used out of the box, that use
[`glam`](https://docs.rs/glam/latest/glam/) to represent the geometry.
*/

use crate::mesh::{
    self, Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, VectorAngleAdaptor,
    VectorLengthAdaptor, VectorNormalizeAdaptor,
};

/// Built-in adaptor for meshes that use 32-bit floating point numbers to
/// represetn the geometry of the mesh.
///
/// This uses [`glam`](https://docs.rs/glam/latest/glam/) to represent the
/// geometry.
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

    fn to_f32(val: Self::Scalar) -> f32 {
        val
    }

    fn to_f64(val: Self::Scalar) -> f64 {
        val as f64
    }
}

/// Built-in adaptor for meshes that use 64-bit floating point numbers to
/// represent the geometry of the mesh.
///
/// This uses [`glam`](https://docs.rs/glam/latest/glam/) to represent the
/// geometry.
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

    fn to_f32(val: Self::Scalar) -> f32 {
        val as f32
    }

    fn to_f64(val: Self::Scalar) -> f64 {
        val
    }
}

/// Mesh type that uses 32 bit floating point numbers to represent the
/// geometry. This mesh implements all adaptors and inherits all features of
/// this crates.
///
/// This uses [`glam`](https://docs.rs/glam/latest/glam/) to represent the
/// geometry.
pub type PolyMeshF32 = mesh::PolyMeshT<3, BuiltInAdaptorF32>;

/// Mesh type that uses 64 bit floating point numbers to represent the
/// geometry. This mesh implements all adaptors and inherits all features of
/// this crates.
///
/// This uses [`glam`](https://docs.rs/glam/latest/glam/) to represent the
/// geometry.
pub type PolyMeshF64 = mesh::PolyMeshT<3, BuiltInAdaptorF64>;
