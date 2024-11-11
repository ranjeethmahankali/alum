pub trait Adaptor<const DIM: usize>
where
    Self::Vector: Default + Clone + Copy,
    Self::Scalar: Default + Clone + Copy,
{
    type Vector;
    type Scalar;

    fn vec(coords: [Self::Scalar; DIM]) -> Self::Vector;

    fn zero_vec() -> Self::Vector;

    fn vec_coord(v: &Self::Vector, i: usize) -> Self::Scalar;
}

pub trait VectorNormalizeAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn normalized_vec(v: Self::Vector) -> Self::Vector;
}

pub trait DotProductAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

pub trait VectorAngleAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

pub trait VectorCrossProductAdaptor<const DIM: usize>: Adaptor<DIM> {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector;
}

pub struct DefaultAdaptorF32 {}

impl Adaptor<3> for DefaultAdaptorF32 {
    type Vector = glam::Vec3;
    type Scalar = f32;

    fn vec(coords: [Self::Scalar; 3]) -> Self::Vector {
        glam::vec3(coords[0], coords[1], coords[2])
    }

    fn zero_vec() -> Self::Vector {
        glam::Vec3::splat(0.)
    }

    fn vec_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
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

impl VectorAngleAdaptor<3> for DefaultAdaptorF32 {
    fn angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl VectorCrossProductAdaptor<3> for DefaultAdaptorF32 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

pub struct DefaultAdaptorF64 {}

impl Adaptor<3> for DefaultAdaptorF64 {
    type Vector = glam::DVec3;
    type Scalar = f64;

    fn vec(coords: [Self::Scalar; 3]) -> Self::Vector {
        glam::dvec3(coords[0], coords[1], coords[2])
    }

    fn zero_vec() -> Self::Vector {
        glam::DVec3::splat(0.)
    }

    fn vec_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
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

impl VectorAngleAdaptor<3> for DefaultAdaptorF64 {
    fn angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.angle_between(b)
    }
}

impl VectorCrossProductAdaptor<3> for DefaultAdaptorF64 {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}
