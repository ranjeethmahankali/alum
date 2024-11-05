use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
    iterator,
    mesh::PolyMeshT,
    property::TPropData,
};
use std::ops::{Add, Div, Mul, Neg, Sub};

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
        Self::Scalar: TScalar + TPropData + Div<Output = Self::Scalar> + PartialOrd,
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

    fn cross(a: Self, b: Self) -> Self
    where
        Self::Scalar: Mul<Output = Self::Scalar> + Sub<Output = Self::Scalar>,
    {
        Self::new(
            a.y() * b.z() - a.z() * b.y(),
            a.z() * b.x() - a.x() * b.z(),
            a.x() * b.y() - a.y() * b.x(),
        )
    }

    fn dot(a: Self, b: Self) -> Self::Scalar
    where
        Self::Scalar: Mul<Output = Self::Scalar> + Add<Output = Self::Scalar>,
    {
        a.x() * b.x() + a.y() * b.y() + a.z() * b.z()
    }

    fn angle(a: Self, b: Self) -> Self::Scalar;
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

    fn normalized(&self) -> Self {
        glam::Vec3::normalize(self.clone())
    }

    fn cross(a: Self, b: Self) -> Self {
        a.cross(b)
    }

    fn dot(a: Self, b: Self) -> Self::Scalar {
        a.dot(b)
    }

    fn angle(a: Self, b: Self) -> Self::Scalar {
        Self::angle_between(a, b)
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT::Scalar: TScalar
        + Mul<Output = VecT::Scalar>
        + Add<Output = VecT::Scalar>
        + TPropData
        + Div<Output = VecT::Scalar>
        + PartialOrd,
    VecT: TVec3 + Add<Output = VecT> + Sub<Output = VecT>,
{
    pub fn try_calc_face_normal(&self, f: FH) -> Result<VecT, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_normal(f, &points))
    }

    pub fn calc_face_normal(&self, f: FH, points: &[VecT]) -> VecT {
        // Use newell's method to compute the normal.
        let (nverts, x, y, z) = {
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
            return VecT::new(
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
            );
        }
        VecT::new(x, y, z).normalized()
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT::Scalar: TScalar
        + Mul<Output = VecT::Scalar>
        + Sub<Output = VecT::Scalar>
        + TPropData
        + Div<Output = VecT::Scalar>
        + PartialOrd,
    VecT: TVec3 + Sub<Output = VecT> + Add<Output = VecT>,
{
    pub fn calc_vertex_normal_accurate(&self, v: VH, points: &[VecT]) -> VecT {
        let topol = self.topology();
        match topol.vertex_halfedge(v) {
            Some(h) => {
                let h2 = topol.ccw_rotated_halfedge(h);
                if h2 == h {
                    // Isolated vertex.
                    return VecT::default();
                }
                // Iterate over adjacent pairs of outgoing halfedges.
                iterator::ccw_rotate_iter(topol, h).zip(iterator::ccw_rotate_iter(topol, h2))
            }
            None => return VecT::default(),
        }
        .fold(
            VecT::new(
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
                VecT::Scalar::from_f64(0.0),
            ),
            |total, (h1, h2)| {
                // Intentionally not normalizing to account for sector area.
                total
                    + VecT::cross(
                        self.calc_halfedge_vector(h1, points),
                        self.calc_halfedge_vector(h2, points),
                    )
            },
        )
        .normalized()
    }

    pub fn try_calc_vertex_normal_accurate(&self, v: VH) -> Result<VecT, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_vertex_normal_accurate(v, &points))
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT::Scalar: TScalar + Add<Output = VecT::Scalar>,
    VecT: TVec3 + Add<Output = VecT> + Div<VecT::Scalar, Output = VecT>,
{
    pub fn try_calc_face_centroid(&self, f: FH) -> Result<VecT, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_centroid(f, &points))
    }

    pub fn calc_face_centroid(&self, f: FH, points: &[VecT]) -> VecT {
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
        total / denom
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Add<Output = VecT>,
    VecT::Scalar: TScalar + TPropData + Div<Output = VecT::Scalar> + PartialOrd,
{
    pub fn calc_vertex_normal_fast(&self, v: VH, fnormals: &[VecT]) -> VecT {
        iterator::vf_ccw_iter(self.topology(), v)
            .fold(
                VecT::new(
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                    VecT::Scalar::from_f64(0.0),
                ),
                |total, f| total + fnormals[f.index() as usize],
            )
            .normalized()
    }

    pub fn try_calc_vertex_normal_fast(&self, v: VH) -> Result<VecT, Error> {
        match self.face_normals() {
            Some(fnormals) => {
                let fnormals = fnormals.try_borrow()?;
                Ok(self.calc_vertex_normal_fast(v, &fnormals))
            }
            None => Err(Error::FaceNormalsNotAvailable),
        }
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Sub<Output = VecT>,
{
    pub fn try_calc_halfedge_vector(&self, h: HH) -> Result<VecT, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_halfedge_vector(h, &points))
    }

    pub fn calc_halfedge_vector(&self, h: HH, points: &[VecT]) -> VecT {
        points[self.to_vertex(h).index() as usize] - points[self.from_vertex(h).index() as usize]
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Sub<Output = VecT>,
    VecT::Scalar: TScalar
        + Mul<Output = VecT::Scalar>
        + Sub<Output = VecT::Scalar>
        + Add<Output = VecT::Scalar>,
{
    pub fn calc_sector_normal(&self, h: HH, points: &[VecT]) -> VecT {
        VecT::cross(
            self.calc_halfedge_vector(self.topology().prev_halfedge(h), points),
            self.calc_halfedge_vector(h, points),
        )
    }

    pub fn try_calc_sector_normal(&self, h: HH) -> Result<VecT, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_sector_normal(h, &points))
    }

    pub fn calc_sector_area(&self, h: HH, points: &[VecT]) -> VecT::Scalar {
        self.calc_sector_normal(h, points).length() * VecT::Scalar::from_f64(0.5)
    }

    pub fn try_calc_sector_area(&self, h: HH) -> Result<VecT::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_sector_area(h, &points))
    }

    pub fn calc_face_area(&self, f: FH, points: &[VecT]) -> VecT::Scalar {
        self.topology().triangulated_face_vertices(f).fold(
            VecT::Scalar::from_f64(0.),
            |total, vs| {
                let p0 = points[vs[0].index() as usize];
                total
                    + VecT::cross(
                        points[vs[1].index() as usize] - p0,
                        points[vs[2].index() as usize] - p0,
                    )
                    .length()
                        * VecT::Scalar::from_f64(0.5)
            },
        )
    }

    pub fn try_calc_face_area(&self, f: FH) -> Result<VecT::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_area(f, &points))
    }

    pub fn calc_area(&self, points: &[VecT]) -> VecT::Scalar {
        self.faces().fold(VecT::Scalar::from_f64(0.), |total, f| {
            total + self.calc_face_area(f, &points)
        })
    }

    pub fn try_calc_area(&self) -> Result<VecT::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_area(&points))
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Sub<Output = VecT>,
    VecT::Scalar: TScalar
        + Mul<Output = VecT::Scalar>
        + Add<Output = VecT::Scalar>
        + Div<Output = VecT::Scalar>
        + Sub<Output = VecT::Scalar>,
{
    pub fn calc_volume(&self, points: &[VecT]) -> VecT::Scalar {
        if self
            .halfedges()
            .any(|h| self.topology().is_boundary_halfedge(h))
        {
            // Not closed.
            return VecT::Scalar::from_f64(0.);
        }
        self.triangulated_vertices()
            .fold(VecT::Scalar::from_f64(0.), |total, vs| {
                let (p0, p1, p2) = (
                    points[vs[0].index() as usize],
                    points[vs[1].index() as usize],
                    points[vs[2].index() as usize],
                );
                total + (VecT::dot(p0, VecT::cross(p1 - p0, p2 - p0)) / VecT::Scalar::from_f64(6.))
            })
    }

    pub fn try_calc_volume(&self) -> Result<VecT::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_volume(&points))
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Sub<Output = VecT>,
    VecT::Scalar: TScalar
        + Mul<Output = VecT::Scalar>
        + Sub<Output = VecT::Scalar>
        + Add<Output = VecT::Scalar>
        + Neg<Output = VecT::Scalar>
        + PartialOrd,
{
    fn aligned_angle(norm0: VecT, norm1: VecT, hvec: VecT) -> VecT::Scalar {
        if VecT::dot(VecT::cross(norm0, norm1), hvec) >= VecT::Scalar::from_f64(0.) {
            VecT::angle(norm0, norm1)
        } else {
            -VecT::angle(norm0, norm1)
        }
    }

    pub fn calc_dihedral_angle(&self, e: EH, points: &[VecT]) -> VecT::Scalar {
        if self.is_boundary_edge(e) {
            return VecT::Scalar::from_f64(0.);
        }
        let h0 = self.edge_halfedge(e, false);
        let h1 = self.edge_halfedge(e, true);
        Self::aligned_angle(
            self.calc_sector_normal(h0, points),
            self.calc_sector_normal(h1, points),
            self.calc_halfedge_vector(h0, points),
        )
    }

    pub fn try_calc_dihedral_angle(&self, e: EH) -> Result<VecT::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_dihedral_angle(e, &points))
    }

    pub fn calc_dihedral_angle_fast(
        &self,
        e: EH,
        points: &[VecT],
        face_normals: &[VecT],
    ) -> VecT::Scalar {
        let h0 = self.edge_halfedge(e, false);
        let h1 = self.edge_halfedge(e, true);
        let f0 = self.halfedge_face(h0);
        let f1 = self.halfedge_face(h1);
        match (f0, f1) {
            (None, None) | (None, Some(_)) | (Some(_), None) => VecT::Scalar::from_f64(0.),
            (Some(f0), Some(f1)) => Self::aligned_angle(
                face_normals[f0.index() as usize],
                face_normals[f1.index() as usize],
                self.calc_halfedge_vector(h0, points),
            ),
        }
    }

    pub fn try_calc_dihedral_angle_fast(&self, e: EH) -> Result<VecT::Scalar, Error> {
        let fnormals = self.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
        let fnormals = fnormals.try_borrow()?;
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_dihedral_angle_fast(e, &points, &fnormals))
    }
}

impl<VecT> PolyMeshT<VecT>
where
    VecT: TVec3 + Sub<Output = VecT>,
{
    pub fn calc_sector_angle(
        &self,
        h: HH,
        points: &[VecT],
        _face_normals: &[VecT],
    ) -> VecT::Scalar {
        let v0 = self.calc_halfedge_vector(h, points);
        let v1 = self.calc_halfedge_vector(self.topology().opposite_halfedge(h), points);
        let _angle = VecT::angle(v0, v1);
        todo!()
    }
}

#[cfg(test)]
mod test {
    use crate::{error::Error, mesh::PolyMeshF32};

    #[test]
    fn t_box_face_normals() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(
            qbox.faces()
                .map(|f| {
                    qbox.try_calc_face_normal(f)
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
                    qbox.try_calc_face_centroid(f)
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
                    qbox.try_calc_halfedge_vector(h)
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

    #[test]
    fn t_box_accurate_vertex_normals() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(
            qbox.vertices()
                .map(|v| qbox
                    .try_calc_vertex_normal_accurate(v)
                    .expect("Cannot compute the correct vertex normal"))
                .collect::<Vec<_>>(),
            &[
                glam::vec3(-0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, 0.57735026, 0.57735026),
                glam::vec3(-0.57735026, 0.57735026, 0.57735026),
            ]
        );
    }

    #[test]
    fn t_box_fast_vertex_normals() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert!(!qbox.has_face_normals());
        // Update and verify face normals.
        let fnormals = qbox
            .update_face_normals()
            .expect("Unable to update face normals");
        assert!(qbox.has_face_normals());
        let fnormals = fnormals.try_borrow().expect("Cannot borrow face normals");
        let fnormals: &[glam::Vec3] = &fnormals;
        for n in fnormals {
            println!("glam::vec3({:.1}, {:.1}, {:.1}),", n.x(), n.y(), n.z());
        }
        assert_eq!(
            &fnormals,
            &[
                glam::vec3(0.0, 0.0, -1.0),
                glam::vec3(0.0, -1.0, 0.0),
                glam::vec3(1.0, 0.0, 0.0),
                glam::vec3(0.0, 1.0, 0.0),
                glam::vec3(-1.0, 0.0, 0.0),
                glam::vec3(0.0, 0.0, 1.0),
            ]
        );
        let qbox = qbox; // Immutable.
                         // Compute and verify vertex normals.
        assert_eq!(
            qbox.vertices()
                .map(|v| qbox.calc_vertex_normal_fast(v, fnormals))
                .collect::<Vec<_>>(),
            &[
                glam::vec3(-0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, 0.57735026, 0.57735026),
                glam::vec3(-0.57735026, 0.57735026, 0.57735026),
            ]
        );
    }

    #[test]
    fn t_box_sector_area() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        let points = qbox.points();
        let points = points
            .try_borrow()
            .expect("Cannot borrow the points property");
        for h in qbox.halfedges() {
            assert_eq!(qbox.calc_sector_area(h, &points), 0.5);
        }
    }

    #[test]
    fn t_box_face_areas() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        let points = qbox.points();
        let points = points.try_borrow().expect("Cannot borrow points");
        for f in qbox.faces() {
            assert_eq!(qbox.calc_face_area(f, &points), 1.);
        }
    }

    #[test]
    fn t_box_total_area() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(qbox.try_calc_area().expect("Unable to calculate area"), 6.);
    }

    #[test]
    fn t_box_volume() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert_eq!(qbox.try_calc_volume().expect("Cannot compute volume"), 1.);
    }

    #[test]
    fn t_box_update_vertex_normals_fast() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert!(!qbox.has_vertex_normals());
        let vnormals = qbox
            .update_vertex_normals_fast()
            .expect("Cannot update vertex normals");
        let vnormals = vnormals.try_borrow().expect("Cannot borrow vertex normals");
        let vnormals: &[glam::Vec3] = &vnormals;
        assert!(qbox.has_vertex_normals());
        assert_eq!(
            vnormals,
            &[
                glam::vec3(-0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, 0.57735026, 0.57735026),
                glam::vec3(-0.57735026, 0.57735026, 0.57735026),
            ]
        );
    }

    #[test]
    fn t_box_update_vertex_normals_accurate() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        assert!(!qbox.has_vertex_normals());
        let vnormals = qbox
            .update_vertex_normals_accurate()
            .expect("Cannot update vertex normals");
        let vnormals = vnormals.try_borrow().expect("Cannot borrow vertex normals");
        let vnormals: &[glam::Vec3] = &vnormals;
        assert!(qbox.has_vertex_normals());
        assert_eq!(
            vnormals,
            &[
                glam::vec3(-0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, -0.57735026, -0.57735026),
                glam::vec3(0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, 0.57735026, -0.57735026),
                glam::vec3(-0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, -0.57735026, 0.57735026),
                glam::vec3(0.57735026, 0.57735026, 0.57735026),
                glam::vec3(-0.57735026, 0.57735026, 0.57735026),
            ]
        );
    }

    #[test]
    fn t_box_sector_normals() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        let points = qbox.points();
        let points = points.try_borrow().expect("Cannot borrow points");
        let points: &[glam::Vec3] = &points;
        for f in qbox.faces() {
            let fnorm = qbox.calc_face_normal(f, points);
            for h in qbox.fh_ccw_iter(f) {
                assert_eq!(qbox.calc_sector_normal(h, points), fnorm);
            }
        }
    }

    #[test]
    fn t_box_fast_dihedral_angles() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        for e in qbox.edges() {
            let out = qbox.try_calc_dihedral_angle_fast(e);
            assert!(matches!(out, Err(Error::FaceNormalsNotAvailable)));
        }
        qbox.update_face_normals()
            .expect("Cannot update face normals");
        let qbox = qbox;
        for e in qbox.edges() {
            assert_eq!(
                1.5707963,
                qbox.try_calc_dihedral_angle_fast(e)
                    .expect("Cannot compute dihedral angle")
            );
        }
    }

    #[test]
    fn t_box_dihedral_angles() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0., 0., 0.), glam::vec3(1., 1., 1.))
            .expect("Cannot create a box primitive");
        for e in qbox.edges() {
            assert_eq!(
                1.5707963,
                qbox.try_calc_dihedral_angle(e)
                    .expect("Cannot compute dihedral angle")
            );
        }
    }
}
