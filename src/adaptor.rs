pub trait Adaptor<const DIM: usize>
where
    Self::Vector: Clone + Copy,
    Self::Scalar: Clone + Copy,
{
    type Vector;
    type Scalar;

    fn vector(coords: [Self::Scalar; DIM]) -> Self::Vector;

    fn zero_vector() -> Self::Vector;

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar;
}

pub trait VectorLengthAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn vector_length(v: Self::Vector) -> Self::Scalar;
}

pub trait VectorNormalizeAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn normalized_vec(v: Self::Vector) -> Self::Vector;
}

pub trait DotProductAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

pub trait VectorAngleAdaptor: Adaptor<3> {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

pub trait CrossProductAdaptor: Adaptor<3> {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector;
}

pub trait FloatScalarAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn scalarf32(val: f32) -> Self::Scalar;

    fn scalarf64(val: f64) -> Self::Scalar;
}

pub struct DefaultAdaptorF32 {}

impl Adaptor<3> for DefaultAdaptorF32 {
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

impl VectorLengthAdaptor<3> for DefaultAdaptorF32 {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        glam::Vec3::length(v)
    }
}

impl VectorNormalizeAdaptor<3> for DefaultAdaptorF32 {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for DefaultAdaptorF32 {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(b)
    }
}

impl VectorAngleAdaptor for DefaultAdaptorF32 {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl CrossProductAdaptor for DefaultAdaptorF32 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

impl FloatScalarAdaptor<3> for DefaultAdaptorF32 {
    fn scalarf32(val: f32) -> Self::Scalar {
        val
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val as f32
    }
}

pub struct DefaultAdaptorF64 {}

impl Adaptor<3> for DefaultAdaptorF64 {
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

impl VectorLengthAdaptor<3> for DefaultAdaptorF64 {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        glam::DVec3::length(v)
    }
}

impl VectorNormalizeAdaptor<3> for DefaultAdaptorF64 {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

impl DotProductAdaptor<3> for DefaultAdaptorF64 {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(b)
    }
}

impl VectorAngleAdaptor for DefaultAdaptorF64 {
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl CrossProductAdaptor for DefaultAdaptorF64 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

impl FloatScalarAdaptor<3> for DefaultAdaptorF64 {
    fn scalarf32(val: f32) -> Self::Scalar {
        val as f64
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val
    }
}
