use super::Decimater;
use crate::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, Error, FloatScalarAdaptor, Handle,
    HasIterators, HasTopology, PolyMeshT, VectorLengthAdaptor, VectorNormalizeAdaptor, HH, VH,
};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

// Credit to:
// https://github.com/Philip-Trettner/probabilistic-quadrics/blob/master/probabilistic-quadrics.hh
struct Quadric<A>
where
    A: Adaptor<3>,
{
    a00: A::Scalar,
    a01: A::Scalar,
    a02: A::Scalar,
    a11: A::Scalar,
    a12: A::Scalar,
    a22: A::Scalar,
    b0: A::Scalar,
    b1: A::Scalar,
    b2: A::Scalar,
    c: A::Scalar,
}

impl<A> Default for Quadric<A>
where
    A: FloatScalarAdaptor<3>,
{
    fn default() -> Self {
        Self {
            a00: A::scalarf64(0.0),
            a01: A::scalarf64(0.0),
            a02: A::scalarf64(0.0),
            a11: A::scalarf64(0.0),
            a12: A::scalarf64(0.0),
            a22: A::scalarf64(0.0),
            b0: A::scalarf64(0.0),
            b1: A::scalarf64(0.0),
            b2: A::scalarf64(0.0),
            c: A::scalarf64(0.0),
        }
    }
}

impl<A> Clone for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: Clone,
{
    fn clone(&self) -> Self {
        Self {
            a00: self.a00.clone(),
            a01: self.a01.clone(),
            a02: self.a02.clone(),
            a11: self.a11.clone(),
            a12: self.a12.clone(),
            a22: self.a22.clone(),
            b0: self.b0.clone(),
            b1: self.b1.clone(),
            b2: self.b2.clone(),
            c: self.c.clone(),
        }
    }
}

impl<A> Copy for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: Copy,
{
}

impl<A> Quadric<A>
where
    A: FloatScalarAdaptor<3>,
{
    fn plane_quadric(pos: A::Vector, normal: A::Vector) -> Self
    where
        A: DotProductAdaptor<3>,
        A::Scalar: Mul<Output = A::Scalar>,
        A::Vector: Mul<A::Scalar, Output = A::Vector>,
    {
        let (nx, ny, nz) = (
            A::vector_coord(&normal, 0),
            A::vector_coord(&normal, 1),
            A::vector_coord(&normal, 2),
        );
        let dot = A::dot_product(pos, normal);
        let ndv = normal * dot;
        Self {
            a00: nx * nx,
            a01: nx * ny,
            a02: nx * nz,
            a11: ny * ny,
            a12: ny * nz,
            a22: nz * nz,
            b0: A::vector_coord(&ndv, 0),
            b1: A::vector_coord(&ndv, 1),
            b2: A::vector_coord(&ndv, 2),
            c: dot * dot,
        }
    }

    fn minimizer(&self) -> A::Vector
    where
        A::Scalar: Div<Output = A::Scalar>
            + Neg<Output = A::Scalar>
            + Mul<Output = A::Scalar>
            + Add<Output = A::Scalar>,
    {
        let one = A::scalarf64(1.0);
        let a00 = self.a00;
        let a10 = self.a01;
        let a20 = self.a02;
        let mut a11 = self.a12;
        let mut a21 = self.a22;
        let mut a22 = self.a22;
        let mut x0 = self.b0;
        let mut x1 = self.b1;
        let mut x2 = self.b2;
        let d0 = one / a00;
        let l10 = a10 * (-d0);
        let l20 = a20 * (-d0);
        a11 = a10 * l10 + a11;
        a21 = a20 * l10 + a21;
        a22 = a20 * l20 + a22;
        let d1 = one / a11;
        let l21 = a21 * (-d1);
        a22 = a21 * l21 + a22;
        let d2 = one / a22;
        x1 = l10 * x0 + x1;
        x2 = l20 * x0 + x2;
        x2 = l21 * x1 + x2;
        x0 = x0 * d0;
        x1 = x1 * d1;
        x2 = x2 * d2;
        x0 = l20 * x2 + x0;
        x1 = l21 * x2 + x1;
        x0 = l10 * x1 + x0;
        A::vector([x0, x1, x2])
    }

    fn residual(&self, p: A::Vector) -> A::Scalar
    where
        A: DotProductAdaptor<3>,
        A::Scalar: Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Sub<Output = A::Scalar>,
    {
        let p0 = A::vector_coord(&p, 0);
        let p1 = A::vector_coord(&p, 1);
        let p2 = A::vector_coord(&p, 2);
        A::dot_product(
            p,
            A::vector([
                self.a00 * p0 + self.a01 * p1 + self.a02 * p2,
                self.a01 * p0 + self.a11 * p1 + self.a12 * p2,
                self.a02 * p0 + self.a12 * p1 + self.a22 * p2,
            ]),
        ) - A::scalarf64(2.0) * (p0 * self.b0 + p1 * self.b1 + p2 * self.b2)
            + self.c
    }
}

impl<A> AddAssign for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: AddAssign,
{
    fn add_assign(&mut self, rhs: Self) {
        self.a00 += rhs.a00;
        self.a01 += rhs.a01;
        self.a02 += rhs.a02;
        self.a11 += rhs.a11;
        self.a12 += rhs.a12;
        self.a22 += rhs.a22;
        self.b0 += rhs.b0;
        self.b1 += rhs.b1;
        self.b2 += rhs.b2;
        self.c += rhs.c;
    }
}

impl<A> Add for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: Add<Output = A::Scalar>,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            a00: self.a00 + rhs.a00,
            a01: self.a01 + rhs.a01,
            a02: self.a02 + rhs.a02,
            a11: self.a11 + rhs.a11,
            a12: self.a12 + rhs.a12,
            a22: self.a22 + rhs.a22,
            b0: self.b0 + rhs.b0,
            b1: self.b1 + rhs.b1,
            b2: self.b2 + rhs.b2,
            c: self.c + rhs.c,
        }
    }
}

impl<A> SubAssign for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: SubAssign,
{
    fn sub_assign(&mut self, rhs: Self) {
        self.a00 -= rhs.a00;
        self.a01 -= rhs.a01;
        self.a02 -= rhs.a02;
        self.a11 -= rhs.a11;
        self.a12 -= rhs.a12;
        self.a22 -= rhs.a22;
        self.b0 -= rhs.b0;
        self.b1 -= rhs.b1;
        self.b2 -= rhs.b2;
        self.c -= rhs.c;
    }
}

impl<A> Neg for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: Neg<Output = A::Scalar>,
{
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self {
            a00: -self.a00,
            a01: -self.a01,
            a02: -self.a02,
            a11: -self.a11,
            a12: -self.a12,
            a22: -self.a22,
            b0: -self.b0,
            b1: -self.b1,
            b2: -self.b2,
            c: -self.c,
        }
    }
}

pub struct QuadricDecimater<A>
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Add<Output = A::Vector> + Sub<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    quadrics: Vec<Quadric<A>>,
}

impl<A> QuadricDecimater<A>
where
    A: CrossProductAdaptor
        + DotProductAdaptor<3>
        + VectorNormalizeAdaptor<3>
        + FloatScalarAdaptor<3>
        + VectorLengthAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector: Add<Output = A::Vector>
        + Sub<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    pub fn new(mesh: &PolyMeshT<3, A>) -> Result<Self, Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        let fnormals = mesh.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
        let fnormals = fnormals.try_borrow()?;
        Ok(Self {
            quadrics: mesh
                .vertices()
                .map(|v| {
                    let pos = points[v.index() as usize];
                    mesh.vf_ccw_iter(v)
                        .fold(Quadric::<A>::default(), |quadric, f| {
                            quadric + Quadric::<A>::plane_quadric(pos, fnormals[f.index() as usize])
                        })
                })
                .collect(),
        })
    }
}

impl<A> Decimater<PolyMeshT<3, A>> for QuadricDecimater<A>
where
    A: CrossProductAdaptor
        + VectorNormalizeAdaptor<3>
        + FloatScalarAdaptor<3>
        + DotProductAdaptor<3>,
    A::Scalar: Mul<Output = A::Scalar>
        + Add<Output = A::Scalar>
        + Div<Output = A::Scalar>
        + PartialOrd
        + Sub<Output = A::Scalar>
        + Clone
        + Copy
        + Neg<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Sub<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    type Cost = A::Scalar;

    fn collapse_cost(&self, mesh: &PolyMeshT<3, A>, h: HH) -> Option<A::Scalar> {
        let q = self.quadrics[h.tail(mesh).index() as usize]
            + self.quadrics[h.head(mesh).index() as usize];
        Some(q.residual(q.minimizer()))
    }

    fn before_collapse(&mut self, _mesh: &PolyMeshT<3, A>, _h: HH) -> Result<(), Error> {
        Ok(())
    }

    fn after_collapse(&mut self, mesh: &PolyMeshT<3, A>, v: VH) -> Result<(), Error> {
        let mut fnormals = mesh.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
        let mut fnormals = fnormals.try_borrow_mut()?;
        let points = mesh.points();
        let points = points.try_borrow()?;
        // Update face normals.
        for f in mesh.vf_ccw_iter(v) {
            fnormals[f.index() as usize] = mesh.calc_face_normal(f, &points);
        }
        // Update quadrics. Use tags to avoid double assignment.
        let mut vstatus = mesh.vertex_status_prop();
        let mut vstatus = vstatus.try_borrow_mut()?;
        for f in mesh.vf_ccw_iter(v) {
            for v in mesh.fv_ccw_iter(f) {
                let vi = v.index() as usize;
                let vs = &mut vstatus[vi];
                if vs.tagged() {
                    continue;
                }
                vs.set_tagged(true);
                let pos = points[vi];
                self.quadrics[vi] =
                    mesh.vf_ccw_iter(v)
                        .fold(Quadric::<A>::default(), |quadric, f| {
                            quadric + Quadric::<A>::plane_quadric(pos, fnormals[f.index() as usize])
                        })
            }
        }
        Ok(())
    }
}
