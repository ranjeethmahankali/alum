use crate::{
    element::{Handle, FH, VH},
    error::Error,
    iterator,
    mesh::PolyMeshT,
    property::TPropData,
};
use std::ops::{Add, Div, Mul, Sub};

pub trait TVec3: TPropData {
    type Scalar;

    fn new(x: Self::Scalar, y: Self::Scalar, z: Self::Scalar) -> Self;

    fn x(&self) -> Self::Scalar;
    fn y(&self) -> Self::Scalar;
    fn z(&self) -> Self::Scalar;

    fn length(&self) -> Self::Scalar;

    fn normalized(&self) -> Self
    where
        Self::Scalar: Div<Output = Self::Scalar> + TPropData + From<f64> + PartialOrd,
    {
        let norm = self.length();
        if norm > Self::Scalar::from(0.0) {
            Self::new(self.x() / norm, self.y() / norm, self.z() / norm)
        } else {
            Self::new(0.0.into(), 0.0.into(), 0.0.into())
        }
    }
}

impl<VecT: TVec3> PolyMeshT<VecT>
where
    VecT::Scalar: From<f64>
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
            let points: &Vec<VecT> = &points;
            iterator::fh_ccw_iter(self.topology(), f).fold(
                (0usize, 0.0.into(), 0.0.into(), 0.0.into()),
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
            return Ok(VecT::new(0.0.into(), 0.0.into(), 0.0.into()));
        }
        Ok(VecT::new(x, y, z).normalized())
    }

    pub fn calc_face_centroid(&self, f: FH) -> Result<VecT, Error> {
        let points = self.points().try_borrow()?;
        let points: &Vec<VecT> = &points;
        let (denom, total) = iterator::fv_ccw_iter(self.topology(), f).fold(
            (0.0.into(), VecT::new(0.0.into(), 0.0.into(), 0.0.into())),
            |(denom, total): (VecT::Scalar, VecT), v: VH| {
                (denom + 1.0.into(), total + points[v.index() as usize])
            },
        );
        Ok(total / denom)
    }
}

#[cfg(test)]
mod test {

    #[test]
    fn t_dodecahedron_face_normals() {
        todo!()
    }

    #[test]
    fn t_dodecahedron_vertex_normals() {
        todo!()
    }

    #[test]
    fn t_box_centroid() {
        todo!()
    }
}
