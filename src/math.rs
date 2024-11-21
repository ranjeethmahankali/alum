use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
    iterator::HasIterators,
    mesh::{
        Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, PolyMeshT,
        VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
    },
    property::{FProperty, VProperty},
    HasTopology, VPropRef,
};
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3> + FloatScalarAdaptor<3>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Add<Output = A::Vector> + Sub<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    /// This is similar to [`Self::calc_face_normal`], except this function
    /// attempts to borrow the necessary properties and return an error if
    /// borrowing fails.
    pub fn try_calc_face_normal(&self, f: FH) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_normal(f, &points))
    }

    /// Compute the face normal using Newell's method. The `points` must
    /// represent the positions of the vertices.
    ///
    /// Calling this function in a hot loop with borrowed points property, can
    /// be faster than [`Self::try_calc_face_normal`] because it avoids repeated
    /// borrows.
    pub fn calc_face_normal(&self, f: FH, points: &VPropRef<A::Vector>) -> A::Vector {
        // Use newell's method to compute the normal.
        let (nverts, x, y, z) = {
            self.fh_ccw_iter(f).fold(
                (
                    0usize,
                    A::scalarf64(0.0),
                    A::scalarf64(0.0),
                    A::scalarf64(0.0),
                ),
                |(nverts, x, y, z): (usize, A::Scalar, A::Scalar, A::Scalar), h| {
                    let (a, b) = {
                        let pc = points[h.tail(self)];
                        let pn = points[h.head(self)];
                        (pc - pn, pc + pn)
                    };
                    (
                        nverts + 1,
                        x + A::vector_coord(&a, 1) * A::vector_coord(&b, 2),
                        y + A::vector_coord(&a, 2) * A::vector_coord(&b, 0),
                        z + A::vector_coord(&a, 0) * A::vector_coord(&b, 1),
                    )
                },
            )
        };
        if nverts < 3 {
            // Guard against degenerate cases.
            return A::zero_vector();
        }
        A::normalized_vec(A::vector([x, y, z]))
    }

    /// Compute the face normals. If the face normals property is not available,
    /// it is initialized before computing the face normals.
    pub fn update_face_normals(&mut self) -> Result<FProperty<A::Vector>, Error> {
        let mut fprop = self.request_face_normals();
        {
            let mut fnormals = fprop.try_borrow_mut()?;
            let fnormals: &mut [A::Vector] = &mut fnormals;
            let points = self.points();
            let points = points.try_borrow()?;
            for f in self.faces() {
                fnormals[f.index() as usize] = self.calc_face_normal(f, &points);
            }
        }
        Ok(fprop)
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + VectorNormalizeAdaptor<3>,
    A::Scalar: Div<Output = A::Scalar> + PartialOrd,
    A::Vector:
        Sub<Output = A::Vector> + Add<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    /// Compute the vertex normals accurately.
    ///
    /// The vertex normal is computed as the average of the normals of the
    /// sectors around the vertex. The `points` argument must be the positions
    /// of the vertices. Calling this function with borrowed `points` in a hot
    /// loop can be faster than [`Self::try_calc_vertex_normal_accurate`] by
    /// avoiding repeated borrows.
    pub fn calc_vertex_normal_accurate(&self, v: VH, points: &VPropRef<A::Vector>) -> A::Vector {
        A::normalized_vec(
            match v.halfedge(self) {
                Some(h) => {
                    let h2 = h.rotate_ccw(self);
                    if h2 == h {
                        // Isolated vertex.
                        return A::zero_vector();
                    }
                    // Iterate over adjacent pairs of outgoing halfedges.
                    self.ccw_rotate_iter(h).zip(self.ccw_rotate_iter(h2))
                }
                None => return A::zero_vector(),
            }
            .fold(A::zero_vector(), |total, (h1, h2)| {
                // Intentionally not normalizing to account for sector area.
                total
                    + A::cross_product(
                        self.calc_halfedge_vector(h1, points),
                        self.calc_halfedge_vector(h2, points),
                    )
            }),
        )
    }

    /// This is similar to [`Self::calc_vertex_normal_accurate`], except this
    /// function will attempt to borrow the required properties and return an
    /// error if borrowing fails.
    pub fn try_calc_vertex_normal_accurate(&self, v: VH) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_vertex_normal_accurate(v, &points))
    }

    /// Compute accurate vertex normals. This can be slower than the fast
    /// approximation. If the vertex normals property is not available, it is
    /// initialized before computing the vertex normals.
    pub fn update_vertex_normals_accurate(&mut self) -> Result<VProperty<A::Vector>, Error> {
        let mut vprop = self.request_vertex_normals();
        {
            let mut vnormals = vprop.try_borrow_mut()?;
            let vnormals: &mut [A::Vector] = &mut vnormals;
            let points = self.points();
            let points = points.try_borrow()?;
            for v in self.vertices() {
                vnormals[v.index() as usize] = self.calc_vertex_normal_accurate(v, &points);
            }
        }
        Ok(vprop)
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: FloatScalarAdaptor<DIM>,
    A::Scalar: Add<Output = A::Scalar>,
    A::Vector: Add<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    /// Similar to [`Self::calc_face_centroid`], except this function attempts
    /// to borrow the necessary properties and returns an error if the borrowing
    /// fails.
    pub fn try_calc_face_centroid(&self, f: FH) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_centroid(f, &points))
    }

    /// Compute the centroid of a face as the average position of the incident vertices.
    ///
    /// `points` must be the positions of the vertices.
    pub fn calc_face_centroid(&self, f: FH, points: &VPropRef<A::Vector>) -> A::Vector {
        let (denom, total) = self.fv_ccw_iter(f).fold(
            (A::scalarf64(0.0), A::zero_vector()),
            |(denom, total): (A::Scalar, A::Vector), v: VH| {
                (denom + A::scalarf64(1.0), total + points[v])
            },
        );
        total / denom
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM> + VectorNormalizeAdaptor<DIM>,
    A::Vector: Add<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
    A::Scalar: Div<Output = A::Scalar> + PartialOrd,
{
    /// Compute the vertex normal as the average of normals of incident
    /// faces. The normals of the incident faces are read from provided
    /// `fnormals`.
    pub fn calc_vertex_normal_fast(&self, v: VH, fnormals: &[A::Vector]) -> A::Vector {
        A::normalized_vec(self.vf_ccw_iter(v).fold(A::zero_vector(), |total, f| {
            total + fnormals[f.index() as usize]
        }))
    }

    /// Similar to [`Self::calc_vertex_normal_fast`] except this function will
    /// attempt to borrow the face normals. If the borrowing fails, or if the
    /// face normals are not available, an error is returned.
    pub fn try_calc_vertex_normal_fast(&self, v: VH) -> Result<A::Vector, Error> {
        match self.face_normals() {
            Some(fnormals) => {
                let fnormals = fnormals.try_borrow()?;
                Ok(self.calc_vertex_normal_fast(v, &fnormals))
            }
            None => Err(Error::FaceNormalsNotAvailable),
        }
    }

    /// Compute a fast approximation of the vertex normals. If the vertex
    /// normals property is not available, it is initialized before computing
    /// the vertex normals.
    pub fn update_vertex_normals_fast(&mut self) -> Result<VProperty<A::Vector>, Error> {
        let mut vprop = self.request_vertex_normals();
        {
            let fnormals = self.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
            let fnormals = fnormals.try_borrow()?;
            let mut vnormals = vprop.try_borrow_mut()?;
            let vnormals: &mut [A::Vector] = &mut vnormals;
            for v in self.vertices() {
                vnormals[v.index() as usize] = self.calc_vertex_normal_fast(v, &fnormals);
            }
        }
        Ok(vprop)
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
{
    /// Similar to [`Self::calc_halfedge_vector`], except this function will
    /// attempt to borrow the necessary properties and return an error if the
    /// borrowing fails.
    pub fn try_calc_halfedge_vector(&self, h: HH) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_halfedge_vector(h, &points))
    }

    /// Compute the vector spanning from the start of the halfedge to the head
    /// of the halfedge.
    pub fn calc_halfedge_vector(&self, h: HH, points: &VPropRef<A::Vector>) -> A::Vector {
        points[h.head(self)] - points[h.tail(self)]
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: Mul<Output = A::Scalar> + Sub<Output = A::Scalar> + Add<Output = A::Scalar>,
{
    /// Compute the normal of a sector, using the given `points` as the
    /// positions of vertices.
    ///
    /// A sector is the triangular region defined by the given halfedge and it's
    /// previous halfedge. Calling this function in a hot loop with borrowed
    /// `points` can be faster than [`Self::try_calc_sector_normal`] by avoiding
    /// repeated borrows.
    pub fn calc_sector_normal(&self, h: HH, points: &VPropRef<A::Vector>) -> A::Vector {
        A::cross_product(
            self.calc_halfedge_vector(h.prev(&self.topol), points),
            self.calc_halfedge_vector(h, points),
        )
    }

    /// Similar to [`Self::calc_sector_normal`], except this function attempts
    /// to borrow the necessary properties and returns an error if borrowing
    /// fails.
    pub fn try_calc_sector_normal(&self, h: HH) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_sector_normal(h, &points))
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + VectorLengthAdaptor<3> + FloatScalarAdaptor<3>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: Mul<Output = A::Scalar> + Sub<Output = A::Scalar> + Add<Output = A::Scalar>,
{
    /// Compute the area of a sector. `points` must be the positions of the vertices.
    ///
    /// A sector is the triangular region defined by the given halfedge and it's
    /// previous halfedge. Calling this function in a hotloop with borrowed
    /// points can be faster than [`Self::try_calc_sector_area`] by avoiding
    /// repeated borrows.
    pub fn calc_sector_area(&self, h: HH, points: &VPropRef<A::Vector>) -> A::Scalar {
        A::vector_length(self.calc_sector_normal(h, points)) * A::scalarf64(0.5)
    }

    /// Similar to [`Self::calc_sector_area`], except this function attempts to
    /// borrow the necessary properties, and returns an error if the borrowing
    /// fails.
    pub fn try_calc_sector_area(&self, h: HH) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_sector_area(h, &points))
    }

    /// Compute the area of a face. `point` must be the positions of vertices.
    ///
    /// For non-planar polygonal faces, the computed area will be
    /// approximate. This is because the area is computed as the sum of
    /// triangles, present in the default triangulation of the face.
    pub fn calc_face_area(&self, f: FH, points: &VPropRef<A::Vector>) -> A::Scalar {
        self.topol
            .triangulated_face_vertices(f)
            .fold(A::scalarf64(0.0), |total, vs| {
                let p0 = points[vs[0]];
                total
                    + A::vector_length(A::cross_product(points[vs[1]] - p0, points[vs[2]] - p0))
                        * A::scalarf64(0.5)
            })
    }

    /// Similar to [`Self::calc_face_area`], except this function will attempt
    /// to borrow the necessary properties, and return an error if the borrowing
    /// fails.
    pub fn try_calc_face_area(&self, f: FH) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_face_area(f, &points))
    }

    /// Compute the total area of this mesh. `points` must be the positions of vertices.
    ///
    /// Calling this function with borrowed `points` property avoids an internal borrow.
    pub fn calc_area(&self, points: &VPropRef<A::Vector>) -> A::Scalar {
        self.faces().fold(A::scalarf64(0.0), |total, f| {
            total + self.calc_face_area(f, points)
        })
    }

    /// Similar to [`Self::calc_area`], except this function tries to borrow the
    /// required properties, and returns an error when borrowing fails.
    ///
    /// ```rust
    /// use alum::use_glam::PolyMeshF32;
    ///
    /// let mesh = PolyMeshF32::unit_box().expect("Cannot create mesh");
    /// assert_eq!(6.0, mesh.try_calc_area().expect("Cannot compute area"));
    /// ```
    pub fn try_calc_area(&self) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_area(&points))
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: VectorLengthAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
{
    /// Compute the length of a mesh edge. `points` must be the positions of the vertices.
    pub fn calc_edge_length(&self, e: EH, points: &VPropRef<A::Vector>) -> A::Scalar {
        A::vector_length(self.calc_halfedge_vector(e.halfedge(false), points))
    }

    /// Similar to [`Self::calc_edge_length`], except this function tries to
    /// borrow the required properties, and returns an error when borrowing
    /// fails.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
    ///     .expect("Cannot create a box primitive");
    /// // All edges of a unit cube must be of length 1.
    /// for e in qbox.edges() {
    ///     assert_eq!(
    ///         1.0,
    ///         qbox.try_calc_edge_length(e)
    ///             .expect("Cannot compute edge length")
    ///     );
    /// }
    /// ```
    pub fn try_calc_edge_length(&self, e: EH) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_edge_length(e, &points))
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: DotProductAdaptor<DIM>,
    A::Vector: Sub<Output = A::Vector>,
{
    /// Compute the squared length of a mesh edge. `points` must be the
    /// positions of the vertices.
    pub fn calc_edge_length_sqr(&self, e: EH, points: &VPropRef<A::Vector>) -> A::Scalar {
        let ev = self.calc_halfedge_vector(e.halfedge(false), points);
        A::dot_product(ev, ev)
    }

    /// Similar to [`Self::calc_edge_length_sqr`], except this function tries to
    /// borrow the required properties, and returns an error when borrowing
    /// fails.
    ///
    /// ```rust
    /// use alum::{use_glam::PolyMeshF32, HasTopology};
    ///
    /// let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
    ///     .expect("Cannot create a box primitive");
    /// // All edges of a unit cube must be of squared length 1.
    /// for e in qbox.edges() {
    ///     assert_eq!(
    ///         1.0,
    ///         qbox.try_calc_edge_length_sqr(e)
    ///             .expect("Cannot compute edge squared length")
    ///     );
    /// }
    /// ```
    pub fn try_calc_edge_length_sqr(&self, e: EH) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_edge_length_sqr(e, &points))
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + DotProductAdaptor<3> + FloatScalarAdaptor<3>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar:
        Mul<Output = A::Scalar> + Add<Output = A::Scalar> + Div<Output = A::Scalar> + AddAssign,
{
    /// Compute the volume of the mesh. `points` must be the positions of the
    /// vertices.
    ///
    /// Calling this function with borrowed `points` property avoids an internal
    /// borrow of properties.
    pub fn calc_volume(&self, points: &VPropRef<A::Vector>) -> A::Scalar {
        if self.halfedges().any(|h| h.is_boundary(&self.topol)) {
            // Not closed.
            return A::scalarf64(0.0);
        }
        self.triangulated_vertices()
            .fold(A::scalarf64(0.0), |total, vs| {
                let (p0, p1, p2) = (points[vs[0]], points[vs[1]], points[vs[2]]);
                total + (A::dot_product(p0, A::cross_product(p1 - p0, p2 - p0)) / A::scalarf64(6.0))
            })
    }

    /// Similar to [`Self::calc_volume`], except this function attempts to
    /// borrow the required properties, and returns an error if the borrowing
    /// fails.
    ///
    /// ```rust
    /// use alum::use_glam::PolyMeshF32;
    ///
    /// let mesh = PolyMeshF32::unit_box().expect("Cannot create mesh");
    /// assert_eq!(1.0, mesh.try_calc_volume().expect("Cannot compute volume"));
    /// ```
    pub fn try_calc_volume(&self) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_volume(&points))
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + DotProductAdaptor<3> + VectorAngleAdaptor + FloatScalarAdaptor<3>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: Mul<Output = A::Scalar>
        + Sub<Output = A::Scalar>
        + Add<Output = A::Scalar>
        + Neg<Output = A::Scalar>
        + PartialOrd
        + AddAssign,
{
    fn aligned_angle(norm0: A::Vector, norm1: A::Vector, align: A::Vector) -> A::Scalar {
        if A::dot_product(A::cross_product(norm0, norm1), align) >= A::scalarf64(0.0) {
            A::vector_angle(norm0, norm1)
        } else {
            -A::vector_angle(norm0, norm1)
        }
    }

    /// Compute the internal dihedral angle at an edge. `points` must be the
    /// positions of the vertices of this mesh.
    ///
    /// Calling this function with borrowed `points` property avoids internal
    /// borrows, avoids errors, and can be faster when called in a fast loop
    /// compared to [`Self::try_calc_dihedral_angle`].
    ///
    /// The the sector normals of the halfedges are used to compute the dihedral
    /// angle. This can be more accurate than
    /// [`Self::calc_dihedral_angle_fast`].
    pub fn calc_dihedral_angle(&self, e: EH, points: &VPropRef<A::Vector>) -> A::Scalar {
        if e.is_boundary(self) {
            return A::scalarf64(0.0);
        }
        let h0 = e.halfedge(false);
        let h1 = e.halfedge(true);
        Self::aligned_angle(
            self.calc_sector_normal(h0, points),
            self.calc_sector_normal(h1, points),
            self.calc_halfedge_vector(h0, points),
        )
    }

    /// Similar to [`Self::calc_dihedral_angle`], except this function tries to
    /// borrow the required properties, and returns an error if the borrowing
    /// fails.
    pub fn try_calc_dihedral_angle(&self, e: EH) -> Result<A::Scalar, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_dihedral_angle(e, &points))
    }

    /// Compute the dihedral angle at an edge by comparing the normals of the
    /// incident faces.
    ///
    /// `points` are the positions of the mesh. Using borrowed `face_normals`
    /// and `points` properties avoids internal borrows, and can be faster when
    /// computing dihedral angles in a hot loop compared to
    /// [`Self::try_calc_dihedral_angle`].
    pub fn calc_dihedral_angle_fast(
        &self,
        e: EH,
        points: &VPropRef<A::Vector>,
        face_normals: &[A::Vector],
    ) -> A::Scalar {
        let h0 = e.halfedge(false);
        let h1 = e.halfedge(true);
        let f0 = h0.face(self);
        let f1 = h1.face(self);
        match (f0, f1) {
            (None, None) | (None, Some(_)) | (Some(_), None) => A::scalarf64(0.0),
            (Some(f0), Some(f1)) => Self::aligned_angle(
                face_normals[f0.index() as usize],
                face_normals[f1.index() as usize],
                self.calc_halfedge_vector(h0, points),
            ),
        }
    }

    /// Similar to [`Self::calc_dihedral_angle_fast`], except this function
    /// attempts to borrow the required properties, and returns an error if the
    /// borrowing fails.
    pub fn try_calc_dihedral_angle_fast(&self, e: EH) -> Result<A::Scalar, Error> {
        let fnormals = self.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
        let fnormals = fnormals.try_borrow()?;
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_dihedral_angle_fast(e, &points, &fnormals))
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + DotProductAdaptor<3> + VectorAngleAdaptor + FloatScalarAdaptor<3>,
    A::Vector: Sub<Output = A::Vector>,
    A::Scalar: Mul<Output = A::Scalar> + Neg<Output = A::Scalar> + AddAssign + PartialOrd,
{
    /// Compute the sector angle.
    ///
    /// A sector associated with a halfedge is the triangular region defined by
    /// the halfedge and it's previous halfedge. `points` are the positions of
    /// mesh vertices. Calling this function with borrowed `points` and
    /// `face_normals` properties avoids internal borrows and can be faster than
    /// [`Self::try_calc_sector_angle`] when computing sector angles in a hot
    /// loop.
    pub fn calc_sector_angle(
        &self,
        h: HH,
        points: &VPropRef<A::Vector>,
        face_normals: &[A::Vector],
    ) -> A::Scalar {
        let n0 = self.calc_halfedge_vector(h, points);
        let h2 = h.prev(self).opposite();
        let n1 = self.calc_halfedge_vector(h2, points);
        let angle = A::vector_angle(n0, n1);
        if let Some(f) = h.opposite().face(self) {
            if h.is_boundary(self)
                && A::dot_product(A::cross_product(n0, n1), face_normals[f.index() as usize])
                    < A::scalarf64(0.0)
            {
                return -angle;
            }
        }
        angle
    }

    /// Similar to [`Self::calc_sector_angle`], except this function tries to
    /// borrow the required properties and returns an error if the borrowing
    /// fails.
    pub fn try_calc_sector_angle(&self, h: HH) -> Result<A::Scalar, Error> {
        let fnormals = self.face_normals().ok_or(Error::FaceNormalsNotAvailable)?;
        let fnormals = fnormals.try_borrow()?;
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_sector_angle(h, &points, &fnormals))
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM> + FloatScalarAdaptor<DIM>,
    A::Vector: Add<Output = A::Vector> + Div<A::Scalar, Output = A::Vector>,
{
    /// Compute the vertex centroid, i.e. the average position of all the vertices.
    ///
    /// `points` must be the positions of vertices.
    pub fn calc_vertex_centroid(&self, points: &VPropRef<A::Vector>) -> A::Vector {
        points.iter().fold(A::zero_vector(), |total, p| total + *p)
            / A::scalarf64(points.len() as f64)
    }

    /// Similar to [`Self::calc_vertex_centroid`] except this function tries to borrow
    /// `points` property and returns an error if the borrowing fails.
    pub fn try_calc_vertex_centroid(&self) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_vertex_centroid(&points))
    }
}

impl<A> PolyMeshT<3, A>
where
    A: CrossProductAdaptor + VectorLengthAdaptor<3> + FloatScalarAdaptor<3>,
    A::Vector: Sub<Output = A::Vector>
        + Add<Output = A::Vector>
        + Div<A::Scalar, Output = A::Vector>
        + Mul<A::Scalar, Output = A::Vector>,
    A::Scalar: Mul<Output = A::Scalar> + Sub<Output = A::Scalar> + Add<Output = A::Scalar>,
{
    /// Compute the area weighted centroid of the mesh.
    ///
    /// `points` must be the positions of vertices.
    pub fn calc_area_centroid(&self, points: &VPropRef<A::Vector>) -> A::Vector {
        let (total, denom) = self.faces().fold(
            (A::zero_vector(), A::scalarf64(0.0)),
            |(total, denom), f| {
                let p = self.calc_face_centroid(f, points);
                let w = self.calc_face_area(f, points);
                (total + p * w, denom + w)
            },
        );
        total / denom
    }

    /// Similar to [`Self::calc_area_centroid`] except this function tries to
    /// borrow the points property and returns an error if the borrowing fails.
    pub fn try_calc_area_centroid(&self) -> Result<A::Vector, Error> {
        let points = self.points();
        let points = points.try_borrow()?;
        Ok(self.calc_area_centroid(&points))
    }
}

#[cfg(all(test, feature = "use_glam"))]
mod test {
    use core::f32;

    use crate::{
        error::Error, iterator::HasIterators, macros::assert_f32_eq, use_glam::PolyMeshF32,
        HasTopology,
    };

    #[test]
    fn t_box_face_normals() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        assert!(!qbox.has_face_normals());
        // Update and verify face normals.
        let fnormals = qbox
            .update_face_normals()
            .expect("Unable to update face normals");
        assert!(qbox.has_face_normals());
        let fnormals = fnormals.try_borrow().expect("Cannot borrow face normals");
        let fnormals: &[glam::Vec3] = &fnormals;
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        let points = qbox.points();
        let points = points.try_borrow().expect("Cannot borrow points");
        for f in qbox.faces() {
            assert_eq!(qbox.calc_face_area(f, &points), 1.0);
        }
    }

    #[test]
    fn t_box_total_area() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        assert_eq!(qbox.try_calc_area().expect("Unable to calculate area"), 6.0);
    }

    #[test]
    fn t_box_volume() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        assert_eq!(qbox.try_calc_volume().expect("Cannot compute volume"), 1.0);
    }

    #[test]
    fn t_box_update_vertex_normals_fast() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        assert!(!qbox.has_vertex_normals());
        qbox.update_face_normals()
            .expect("Cannot update face normals");
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
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
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
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        let points = qbox.points();
        let points = points.try_borrow().expect("Cannot borrow points");
        for f in qbox.faces() {
            let fnorm = qbox.calc_face_normal(f, &points);
            for h in qbox.fh_ccw_iter(f) {
                assert_eq!(qbox.calc_sector_normal(h, &points), fnorm);
            }
        }
    }

    #[test]
    fn t_box_fast_dihedral_angles() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        for e in qbox.edges() {
            let out = qbox.try_calc_dihedral_angle_fast(e);
            assert!(matches!(out, Err(Error::FaceNormalsNotAvailable)));
        }
        qbox.update_face_normals()
            .expect("Cannot update face normals");
        let qbox = qbox;
        for e in qbox.edges() {
            assert_f32_eq!(
                f32::consts::FRAC_PI_2,
                qbox.try_calc_dihedral_angle_fast(e)
                    .expect("Cannot compute dihedral angle")
            );
        }
    }

    #[test]
    fn t_box_dihedral_angles() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        for e in qbox.edges() {
            assert_f32_eq!(
                f32::consts::FRAC_PI_2,
                qbox.try_calc_dihedral_angle(e)
                    .expect("Cannot compute dihedral angle")
            );
        }
    }

    #[test]
    fn t_box_sector_angles() {
        let mut qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        assert!(matches!(
            qbox.try_calc_sector_angle(0.into()),
            Err(Error::FaceNormalsNotAvailable)
        ));
        qbox.update_face_normals()
            .expect("Cannot update face normals");
        for h in qbox.halfedges() {
            assert_f32_eq!(
                f32::consts::FRAC_PI_2,
                qbox.try_calc_sector_angle(h)
                    .expect("Cannot compute sector angles")
            );
        }
        // Remove faces and test boundary edges.
        let mut qbox = qbox;
        for fi in 0..3 {
            qbox.delete_face(fi.into(), true)
                .expect("Cannot delete a face");
        }
        qbox.garbage_collection().expect("Cannot garbage collect");
        let qbox = qbox;
        let mut convex = 0usize;
        let mut concave = 0usize;
        for h in qbox.halfedges().filter(|h| h.is_boundary(&qbox)) {
            let angle = qbox
                .try_calc_sector_angle(h)
                .expect("Cannot compute sector angle");
            assert_f32_eq!(f32::consts::FRAC_PI_2, angle.abs());
            if angle < 0. {
                concave += 1;
            } else {
                convex += 1;
            };
        }
        assert_eq!(3, convex);
        assert_eq!(3, concave);
    }

    #[test]
    fn t_box_edge_lengths() {
        let qbox = PolyMeshF32::quad_box(glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 1.0, 1.0))
            .expect("Cannot create a box primitive");
        for e in qbox.edges() {
            assert_eq!(
                1.0,
                qbox.try_calc_edge_length(e)
                    .expect("Cannot compute edge length")
            );
        }
    }
}
