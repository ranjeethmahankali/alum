use crate::{
    element::{Handle, FH, HH, VH},
    error::Error,
    iterator,
    mesh::PolyMeshT,
    property::TPropData,
};
use std::ops::{Add, Div, Mul, Sub};

pub trait TScalar {
    fn from_f64(val: f64) -> Self;

    fn from_f32(val: f32) -> Self;
}

impl TScalar for f32 {
    fn from_f64(val: f64) -> Self {
        val as f32
    }

    fn from_f32(val: f32) -> Self {
        val
    }
}

impl TScalar for f64 {
    fn from_f64(val: f64) -> Self {
        val
    }

    fn from_f32(val: f32) -> Self {
        val as f64
    }
}

pub trait TVec3: TPropData {
    type Scalar;

    fn new(x: Self::Scalar, y: Self::Scalar, z: Self::Scalar) -> Self;

    fn x(&self) -> Self::Scalar;
    fn y(&self) -> Self::Scalar;
    fn z(&self) -> Self::Scalar;

    fn length(&self) -> Self::Scalar;

    fn normalized(&self) -> Self
    where
        Self::Scalar: Div<Output = Self::Scalar> + TPropData + TScalar + PartialOrd,
    {
        let norm = self.length();
        if norm > Self::Scalar::from_f64(0.0) {
            Self::new(self.x() / norm, self.y() / norm, self.z() / norm)
        } else {
            Self::new(
                Self::Scalar::from_f64(0.0),
                Self::Scalar::from_f64(0.0),
                Self::Scalar::from_f64(0.0),
            )
        }
    }
}

impl TVec3 for glam::Vec3 {
    type Scalar = f32;

    fn new(x: Self::Scalar, y: Self::Scalar, z: Self::Scalar) -> Self {
        glam::vec3(x, y, z)
    }

    fn x(&self) -> Self::Scalar {
        self[0]
    }

    fn y(&self) -> Self::Scalar {
        self[1]
    }

    fn z(&self) -> Self::Scalar {
        self[2]
    }

    fn length(&self) -> Self::Scalar {
        glam::Vec3::length(self.clone())
    }
}

impl<VecT: TVec3> PolyMeshT<VecT>
where
    VecT::Scalar: TScalar
        + TPropData
        + Mul<Output = VecT::Scalar>
        + Add<Output = VecT::Scalar>
        + PartialOrd
        + Div<Output = VecT::Scalar>,
    VecT: Add<Output = VecT> + Sub<Output = VecT> + Div<VecT::Scalar, Output = VecT>,
{
    pub fn calc_face_normal(&self, f: FH) -> Result<VecT, Error> {
        // Use newell's method to compute the normal.
        let (nverts, x, y, z) = {
            // Borrow points inside the scope to limit their lifetime.
            let points = self.points().try_borrow()?;
            let points: &[VecT] = &points;
            iterator::fh_ccw_iter(self.topology(), f).fold(
                (
                    0usize,
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                ),
                |(nverts, x, y, z): (usize, VecT::Scalar, VecT::Scalar, VecT::Scalar), h| {
                    let (a, b) = {
                        let pc = points[self.from_vertex(h).index() as usize];
                        let pn = points[self.to_vertex(h).index() as usize];
                        (pc - pn, pc + pn)
                    };
                    (
                        nverts + 1,
                        x + a.y() * b.z(),
                        y + a.z() * b.x(),
                        z + a.x() * b.y(),
                    )
                },
            )
        };
        if nverts < 3 {
            // Guard against degenerate cases.
            return Ok(VecT::new(
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
            ));
        }
        Ok(VecT::new(x, y, z).normalized())
    }

    pub fn calc_face_centroid(&self, f: FH) -> Result<VecT, Error> {
        let points = self.points().try_borrow()?;
        let points: &[VecT] = &points;
        let (denom, total) = iterator::fv_ccw_iter(self.topology(), f).fold(
            (
                VecT::Scalar::from_f64(0.0),
                VecT::new(
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                ),
            ),
            |(denom, total): (VecT::Scalar, VecT), v: VH| {
                (
                    denom + VecT::Scalar::from_f64(1.0),
                    total + points[v.index() as usize],
                )
            },
        );
        Ok(total / denom)
    }

    pub fn calc_halfedge_vector(&self, h: HH) -> Result<VecT, Error> {
        let points = self.points().try_borrow()?;
        let points: &[VecT] = &points;
        Ok(points[self.to_vertex(h).index() as usize]
            - points[self.from_vertex(h).index() as usize])
    }
}

#[cfg(test)]
mod test {
    use crate::mesh::PolyMeshF32;

    #[test]
    fn t_box_face_normals() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(
            qbox.faces()
                .map(|f| {
                    qbox.calc_face_normal(f)
                        .expect("Cannot compute face centroid")
                })
                .collect::<Vec<_>>(),
            &[
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(0.0, -1.0, 0.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(0.0, 0.0, 1.0),
            ]
        );
    }

    #[test]
    fn t_box_face_centroids() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(
            qbox.faces()
                .map(|f| {
                    qbox.calc_face_centroid(f)
                        .expect("Cannot compute face centroid")
                })
                .collect::<Vec<_>>(),
            &[
                glam::vec3(0.5, 0.5, 0.0),
                glam::vec3(0.5, 0.0, 0.5),
                glam::vec3(1.0, 0.5, 0.5),
                glam::vec3(0.5, 1.0, 0.5),
                glam::vec3(0.0, 0.5, 0.5),
                glam::vec3(0.5, 0.5, 1.0),
            ]
        );
    }

    #[test]
    fn t_box_edge_vectors() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(
            qbox.halfedges()
                .map(|h| {
                    qbox.calc_halfedge_vector(h)
                        .expect("Cannot compute halfedge vector")
                })
                .collect::<Vec<_>>(),
            &[
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(0.0, -1.0, 0.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(0.0, -1.0, 0.0),
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(0.0, 0.0, 1.0),
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(0.0, 0.0, 1.0),
                glam::vec3(0.0, 0.0, 1.0),
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(0.0, -1.0, 0.0),
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(0.0, 0.0, 1.0),
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(0.0, -1.0, 0.0)
            ]
        );
    }
}
