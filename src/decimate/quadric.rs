use super::Decimater;
use crate::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, Error, FloatScalarAdaptor, Handle,
    HasIterators, HasTopology, PolyMeshT, VectorLengthAdaptor, VectorNormalizeAdaptor, HH, VH,
};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

// Credit to:
// https://github.com/Philip-Trettner/probabilistic-quadrics/blob/master/probabilistic-quadrics.hh
pub struct Quadric<A>
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

impl<A> PartialEq for Quadric<A>
where
    A: Adaptor<3>,
    A::Scalar: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        self.a00 == other.a00
            && self.a01 == other.a01
            && self.a02 == other.a02
            && self.a11 == other.a11
            && self.a12 == other.a12
            && self.a22 == other.a22
            && self.b0 == other.b0
            && self.b1 == other.b1
            && self.b2 == other.b2
            && self.c == other.c
    }
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
        *self
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
    pub fn plane_quadric(pos: A::Vector, normal: A::Vector) -> Self
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

    pub fn triangle_quadric(p: A::Vector, q: A::Vector, r: A::Vector) -> Self
    where
        A: CrossProductAdaptor + DotProductAdaptor<3>,
        A::Vector: Add<Output = A::Vector>,
        A::Scalar: Mul<Output = A::Scalar>,
    {
        let pxq = A::cross_product(p, q);
        let qxr = A::cross_product(q, r);
        let rxp = A::cross_product(r, p);
        let xsum = pxq + qxr + rxp;
        let det = A::dot_product(pxq, r);
        {
            let x = A::vector_coord(&xsum, 0);
            let y = A::vector_coord(&xsum, 1);
            let z = A::vector_coord(&xsum, 2);
            Self {
                a00: x * x,
                a01: x * y,
                a02: x * z,
                a11: y * y,
                a12: y * z,
                a22: z * z,
                b0: x * det,
                b1: y * det,
                b2: z * det,
                c: det * det,
            }
        }
    }

    fn cross_product_squared_transpose(v: A::Vector) -> [A::Scalar; 6]
    where
        A::Scalar: Add<Output = A::Scalar> + Neg<Output = A::Scalar> + Mul<Output = A::Scalar>,
    {
        let a = A::vector_coord(&v, 0);
        let b = A::vector_coord(&v, 1);
        let c = A::vector_coord(&v, 2);
        let a2 = a * a;
        let b2 = b * b;
        let c2 = c * c;
        [b2 + c2, -a * b, -a * c, a2 + c2, -b * c, a2 + b2]
    }

    pub fn probabilistic_triangle_quadric(
        p: A::Vector,
        q: A::Vector,
        r: A::Vector,
        stddev: A::Scalar,
    ) -> Self
    where
        A: CrossProductAdaptor + DotProductAdaptor<3>,
        A::Vector:
            Sub<Output = A::Vector> + Add<Output = A::Vector> + Mul<A::Scalar, Output = A::Vector>,
        A::Scalar: Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Neg<Output = A::Scalar>,
    {
        let sigma = stddev * stddev;
        let pxq = A::cross_product(p, q);
        let qxr = A::cross_product(q, r);
        let rxp = A::cross_product(r, p);
        let det_pqr = A::dot_product(pxq, r);
        let cross_pqr = pxq + qxr + rxp;
        let pmq = p - q;
        let qmr = q - r;
        let rmp = r - p;
        let ss = sigma * sigma;
        let ss6 = ss * A::scalarf64(6.0);
        let ss2 = ss * A::scalarf64(2.0);
        let amat = {
            let cout = {
                let x = A::vector_coord(&cross_pqr, 0);
                let y = A::vector_coord(&cross_pqr, 1);
                let z = A::vector_coord(&cross_pqr, 2);
                [x * x, x * y, x * z, y * y, y * z, z * z]
            };
            let pmat = Self::cross_product_squared_transpose(pmq);
            let qmat = Self::cross_product_squared_transpose(qmr);
            let rmat = Self::cross_product_squared_transpose(rmp);
            let mut amat = [A::scalarf64(0.); 6];
            for i in 0..6 {
                amat[i] = cout[i] + (pmat[i] + qmat[i] + rmat[i]) * sigma;
            }
            amat[0] = amat[0] + ss6;
            amat[3] = amat[3] + ss6;
            amat[5] = amat[5] + ss6;
            amat
        };
        let b = (cross_pqr * det_pqr)
            - (A::cross_product(pmq, pxq)
                + A::cross_product(qmr, qxr)
                + A::cross_product(rmp, rxp))
                * sigma
            + (p + q + r) * ss2;
        let c = (det_pqr * det_pqr)
            + (A::dot_product(pxq, pxq) + A::dot_product(qxr, qxr) + A::dot_product(rxp, rxp))
                * sigma
            + (A::dot_product(p, p) + A::dot_product(q, q) + A::dot_product(r, r)) * ss2
            + ss6 * sigma;
        Self {
            a00: amat[0],
            a01: amat[1],
            a02: amat[2],
            a11: amat[3],
            a12: amat[4],
            a22: amat[5],
            b0: A::vector_coord(&b, 0),
            b1: A::vector_coord(&b, 1),
            b2: A::vector_coord(&b, 2),
            c,
        }
    }

    pub fn minimizer(&self) -> A::Vector
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
        let mut a11 = self.a11;
        let mut a21 = self.a12;
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

    pub fn residual(&self, p: A::Vector) -> A::Scalar
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

/// This [`Decimater`](Decimater<MeshT>) implementation prioritizes collapses
/// that minimize the quadric error.
///
/// The type of quadric this decimater uses is determined by how it was
/// initialized. See [this](QuadricDecimater::new) for more information. This
/// uses [this
/// implementation](https://github.com/Philip-Trettner/probabilistic-quadrics)
/// by Philip Trettner.
pub struct QuadricDecimater<A>
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Add<Output = A::Vector> + Sub<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    quadrics: Vec<Quadric<A>>,
    next: Quadric<A>,
}

/// Type of quadric error to be used by the [`QuadricDecimater`].
pub enum QuadricType {
    /// Triangle quadric.
    Triangle,
    /// Probabilistic triangle quadric.
    ProbabilisticTriangle,
}

impl<A> QuadricDecimater<A>
where
    A: CrossProductAdaptor
        + DotProductAdaptor<3>
        + VectorNormalizeAdaptor<3>
        + FloatScalarAdaptor<3>
        + VectorLengthAdaptor<3>,
    A::Scalar: Mul<Output = A::Scalar>
        + Add<Output = A::Scalar>
        + Div<Output = A::Scalar>
        + PartialOrd
        + Neg<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector>
        + Sub<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
{
    /// Create a new QuadricDecimater.
    ///
    /// `qtype` determines the type of quadric error used to prioritize edge
    /// collapses.
    pub fn new(mesh: &PolyMeshT<3, A>, qtype: QuadricType) -> Result<Self, Error> {
        let points = mesh.points();
        let points = points.try_borrow()?;
        const EDGE_RATIO: f64 = 0.1;
        let fqs: Vec<_> = mesh
            .faces()
            .map(|f| {
                mesh.triangulated_face_vertices(f)
                    .fold(Quadric::<A>::default(), |q, vs| {
                        q + match qtype {
                            QuadricType::Triangle => Quadric::<A>::triangle_quadric(
                                points[vs[0]],
                                points[vs[1]],
                                points[vs[2]],
                            ),
                            QuadricType::ProbabilisticTriangle => {
                                Quadric::<A>::probabilistic_triangle_quadric(
                                    points[vs[0]],
                                    points[vs[1]],
                                    points[vs[2]],
                                    A::scalarf64(EDGE_RATIO)
                                        * mesh.fe_ccw_iter(f).fold(A::scalarf64(0.), |total, e| {
                                            total + mesh.calc_edge_length(e, &points)
                                        })
                                        / A::scalarf64(f.valence(mesh) as f64),
                                )
                            }
                        }
                    })
            })
            .collect();
        Ok(Self {
            quadrics: mesh
                .vertices()
                .map(|v| {
                    mesh.vf_ccw_iter(v)
                        .fold(Quadric::<A>::default(), |quadric, f| {
                            quadric + fqs[f.index() as usize]
                        })
                })
                .collect(),
            next: Default::default(),
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

    /// The cost of the collapse, i.e. the quadric error.
    fn collapse_cost(&self, mesh: &PolyMeshT<3, A>, h: HH) -> Option<A::Scalar> {
        let q = self.quadrics[h.tail(mesh).index() as usize]
            + self.quadrics[h.head(mesh).index() as usize];
        Some(q.residual(q.minimizer()))
    }

    /// Compute the sum of the quadrics from vertices about to be
    /// collapsed. This is used in [`after_collapse`](Self::after_collapse).
    fn before_collapse(&mut self, mesh: &PolyMeshT<3, A>, h: HH) -> Result<(), Error> {
        self.next = self.quadrics[h.tail(mesh).index() as usize]
            + self.quadrics[h.head(mesh).index() as usize];
        Ok(())
    }

    /// The vertex quadric is updated to be the sum of the two vertices before
    /// collapse, and it is moved to a position that minimizes this quadric.
    fn after_collapse(&mut self, mesh: &PolyMeshT<3, A>, v: VH) -> Result<(), Error> {
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut()?;
        points[v] = self.next.minimizer();
        self.quadrics[v.index() as usize] = self.next;
        Ok(())
    }
}

#[cfg(test)]
mod test {
    use crate::{
        decimate::quadric::QuadricType, obj::test::bunny_mesh, HasDecimation, HasIterators,
        HasTopology, QuadricDecimater,
    };

    #[test]
    fn t_bunny_quadrics() {
        let mut mesh = bunny_mesh();
        let mut decimater = QuadricDecimater::new(&mesh, QuadricType::Triangle).unwrap();
        let (nv, ne, nf) = (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces());
        let nloops_before = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        mesh.decimate(&mut decimater, 500).unwrap();
        mesh.garbage_collection()
            .expect("Garbage collection failed");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(nv - 500, mesh.num_vertices());
        assert!(mesh.num_edges() < ne);
        assert!(mesh.num_faces() < nf);
        let nloops_after = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        assert_eq!(nloops_before, nloops_after);
    }

    #[test]
    fn t_bunny_probabilistic_quadrics() {
        let mut mesh = bunny_mesh();
        let mut decimater =
            QuadricDecimater::new(&mesh, QuadricType::ProbabilisticTriangle).unwrap();
        let (nv, ne, nf) = (mesh.num_vertices(), mesh.num_edges(), mesh.num_faces());
        let nloops_before = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        mesh.decimate(&mut decimater, 500).unwrap();
        mesh.garbage_collection()
            .expect("Garbage collection failed");
        mesh.check_topology().expect("Topological errors found");
        assert_eq!(nv - 500, mesh.num_vertices());
        assert!(mesh.num_edges() < ne);
        assert!(mesh.num_faces() < nf);
        let nloops_after = {
            let mut visited = mesh.create_halfedge_prop(false);
            let mut visited = visited.try_borrow_mut().unwrap();
            mesh.halfedges().fold(0usize, |count, h| {
                if !h.is_boundary(&mesh) || visited[h] {
                    count
                } else {
                    for h2 in mesh.loop_ccw_iter(h) {
                        visited[h2] = true;
                    }
                    count + 1
                }
            })
        };
        assert_eq!(nloops_before, nloops_after);
    }
}
