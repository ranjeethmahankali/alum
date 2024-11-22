use crate::{
    element::{VRange, EH, FH, VH},
    error::Error,
    property::{FProperty, VProperty},
    topol::{HasTopology, TopolCache, Topology},
};

/// Trait for an adaptor that tells this crate how to work with user specified
/// geometric types.
///
/// Implementing an adaptor allows the use of this crate with any geometric
/// type. An implementation of this trait must specify an associated type
/// `Vector` to represent the positions of vertices, and an associated type
/// `Scalar` to represent the coordinates of `Vector` etc. In addition to these
/// associated types, an implementation is expected to provide basic functions
/// to work with `Vector` instances.
pub trait Adaptor<const DIM: usize>
where
    Self::Vector: Clone + Copy,
    Self::Scalar: Clone + Copy,
{
    type Vector;
    type Scalar;

    /// Construct a vector from an array containing the coordinates.
    fn vector(coords: [Self::Scalar; DIM]) -> Self::Vector;

    /// Construct a zero vector.
    fn zero_vector() -> Self::Vector;

    /// Get the coordiate at the given index from a vector.
    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// compute the length of its vector type.
///
/// Implementing this trait is required by functions of the polygon mesh such as
/// [`PolyMeshT<DIM, A>::calc_edge_length`], and [`PolyMeshT<DIM,
/// A>::calc_face_area`] etc.
pub trait VectorLengthAdaptor<const DIM: usize>: Adaptor<DIM> {
    /// Compute the length (i.e. magnitude) of a vector.
    fn vector_length(v: Self::Vector) -> Self::Scalar;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// normalize a vector to have a unit length.
///
/// Implementing this trait is required by functions of the polygon mesh such as
/// [`PolyMeshT<DIM, A>::calc_face_normal`], [`PolyMeshT<DIM,
/// A>::calc_vertex_normal_accurate`], and [`PolyMeshT<DIM,
/// A>::calc_vertex_normal_fast`] etc.
pub trait VectorNormalizeAdaptor<const DIM: usize>: Adaptor<DIM> {
    /// Normalize a vector and return a unit length vector.
    fn normalized_vec(v: Self::Vector) -> Self::Vector;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// compute the dot product of two vectors.
///
/// Implementing this trait is required by functions of the polygon mesh such as
/// [`PolyMeshT<DIM, A>::calc_volume`], [`PolyMeshT<DIM,
/// A>::calc_sector_angle`], and [`PolyMeshT<DIM, A>::calc_dihedral_angle`] etc.
pub trait DotProductAdaptor<const DIM: usize>: Adaptor<DIM> {
    /// Compute the dot product between two vectors.
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// compute angle between its vectors.
///
/// Implementing this trait is required by functions of the polygon mesh such as
/// [`PolyMeshT<DIM, A>::calc_sector_angle`], and [`PolyMeshT<DIM,
/// A>::calc_dihedral_angle`] etc.
pub trait VectorAngleAdaptor: Adaptor<3> {
    /// Compute the angle between two vectors.
    fn vector_angle(a: Self::Vector, b: Self::Vector) -> Self::Scalar;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// compute the cross product of its vectors.
///
/// This is only meant to be used with 3-d vectors, hence the generic `DIM`
/// parameter is constrainted to 3. Implementing this trait is required by
/// functions of the polygon mesh such as [`PolyMeshT<DIM,
/// A>::calc_face_normal`], and [`PolyMeshT<DIM,
/// A>::calc_vertex_normal_accurate`] etc.
pub trait CrossProductAdaptor: Adaptor<3> {
    /// Compute the cross product of two vectors.
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector;
}

/// An adaptor can optionally implement this trait to tell this crate how to
/// instantiate its scalar type from [`f32`] and [`f64`].
///
/// Implementing this trait is required by functions of the polygon mesh such as
/// [`PolyMeshT<DIM, A>::calc_volume`], [`PolyMeshT<DIM, A>::load_obj`], and
/// several others.
pub trait FloatScalarAdaptor<const DIM: usize>: Adaptor<DIM> {
    /// Create a scalar from [`f32`].
    fn scalarf32(val: f32) -> Self::Scalar;

    /// Create a scalar from [`f64`].
    fn scalarf64(val: f64) -> Self::Scalar;
}

/// A polygon mesh.
///
/// `DIM` is the number of spatial dimensions this mesh lives in. This is the
/// same as the number of coordinates the positions of it's vertices are
/// expected to have. `A` must be an implementation of the [`Adaptor`] that
/// tells the mesh type about the geometric types like vectors and scalars and
/// how to work with them.
pub struct PolyMeshT<const DIM: usize, A>
where
    A: Adaptor<DIM>,
    A::Vector: 'static,
    A::Scalar: 'static,
{
    pub(crate) topol: Topology,
    pub(crate) cache: TopolCache,
    points: VProperty<A::Vector>,
    vnormals: Option<VProperty<A::Vector>>,
    fnormals: Option<FProperty<A::Vector>>,
}

impl<const DIM: usize, A> Default for PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Create a new empty mesh.
    fn default() -> Self {
        Self::new()
    }
}

impl<const DIM: usize, A> PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Create a new empty mesh.
    pub fn new() -> Self {
        let mut topol = Topology::new();
        let points = topol.create_vertex_prop(A::zero_vector());
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    /// Create an empty mesh with memory reserved for the given number of
    /// elements.
    ///
    /// The memory is also reserved for all the built-in properties. Normals are
    /// not computed.
    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let mut topol = Topology::with_capacity(nverts, nedges, nfaces);
        let points = VProperty::with_capacity(nverts, &mut topol.vprops, A::zero_vector());
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    /// The position of a vertex.
    pub fn point(&self, vi: VH) -> Result<A::Vector, Error> {
        self.points.get_cloned(vi)
    }

    /// Set the position of a vertex.
    pub fn set_point(&mut self, v: VH, pos: A::Vector) -> Result<(), Error> {
        self.points.set(v, pos)
    }

    /// Get the property corresponding to the positions of the vertices.
    pub fn points(&self) -> VProperty<A::Vector> {
        self.points.clone()
    }

    /// Check if the vertex normals are available.
    pub fn has_vertex_normals(&self) -> bool {
        self.vnormals.is_some()
    }

    /// Get the vertex normals, if available.
    pub fn vertex_normals(&self) -> Option<VProperty<A::Vector>> {
        self.vnormals.clone()
    }

    /// Get the vertex normals property.
    ///
    /// If the vertex normals are not available, a new property for vertex
    /// normals is created and returned. The vertex normals are not
    /// computed. The values are default initialized to zero vectors.
    /// ```rust
    /// use alum::use_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot create a tetrahedron");
    /// let vnormals = mesh.request_vertex_normals();
    /// let vnormals = vnormals.try_borrow().expect("Cannot borrow property");
    /// for v in vnormals.iter() {
    ///     assert_eq!(glam::vec3(0.0, 0.0, 0.0), *v);
    /// }
    /// ```
    pub fn request_vertex_normals(&mut self) -> VProperty<A::Vector> {
        self.vnormals
            .get_or_insert_with(|| self.topol.create_vertex_prop(A::zero_vector()))
            .clone()
    }

    /// Check if the face normals are available.
    pub fn has_face_normals(&self) -> bool {
        self.fnormals.is_some()
    }

    /// Get the face normals, if available.
    pub fn face_normals(&self) -> Option<FProperty<A::Vector>> {
        self.fnormals.clone()
    }

    /// Get the face normals property.
    ///
    /// If the face normals property is not available, a new property for face
    /// normals is created and returned. The face normals are not computed. The
    /// values in the property are default initialized to zero vectors.
    /// ```rust
    /// use alum::use_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot create a tetrahedron");
    /// let fnormals = mesh.request_face_normals();
    /// let fnormals = fnormals.try_borrow().expect("Cannot borrow property");
    /// for v in fnormals.iter() {
    ///     assert_eq!(glam::vec3(0.0, 0.0, 0.0), *v);
    /// }
    /// ```
    pub fn request_face_normals(&mut self) -> FProperty<A::Vector> {
        self.fnormals
            .get_or_insert_with(|| self.topol.create_face_prop(A::zero_vector()))
            .clone()
    }

    /// Add a vertex to this mesh at the given position.
    pub fn add_vertex(&mut self, pos: A::Vector) -> Result<VH, Error> {
        let vi = self.topol.add_vertex()?;
        self.points.set(vi, pos)?;
        Ok(vi)
    }

    /// Add several vertices at once, at the supplied positions.
    ///
    /// If successful, the range of indices of the newly added vertices is
    /// returned.
    /// ```rust
    /// use alum::use_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// let verts = mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// assert_eq!(verts, 0..4);
    /// let verts = [glam::vec3(0.0, 0.0, 1.0), glam::vec3(1.0, 0.0, 1.0),
    ///              glam::vec3(1.0, 1.0, 1.0), glam::vec3(0.0, 1.0, 1.0),];
    /// let verts = mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// assert_eq!(verts, 4..8);
    /// ```
    pub fn add_vertices(&mut self, pos: &[A::Vector]) -> Result<VRange, Error> {
        let vis = self.topol.add_vertices(pos.len())?;
        let mut points = self.points.try_borrow_mut()?;
        for (i, v) in vis.clone().enumerate() {
            points[v] = pos[i];
        }
        Ok(vis)
    }

    /// Add a face to this mesh. The `verts` are expected to be in counter-clockwise order.
    pub fn add_face(&mut self, verts: &[VH]) -> Result<FH, Error> {
        self.topol.add_face(verts, &mut self.cache)
    }

    /// Add a triangular face with the given vertices.
    pub fn add_tri_face(&mut self, v0: VH, v1: VH, v2: VH) -> Result<FH, Error> {
        self.add_face(&[v0, v1, v2])
    }

    /// Add a quadrilateral face with the given vertices.
    pub fn add_quad_face(&mut self, v0: VH, v1: VH, v2: VH, v3: VH) -> Result<FH, Error> {
        self.add_face(&[v0, v1, v2, v3])
    }

    /// Delete a vertex.
    ///
    /// This will automatically delete all incident edges and faces. Vertices in
    /// the neighborhood may become isolated. They are deleted if
    /// `delete_isolated_vertices` is true. The elements are marked as deleted,
    /// but not removed from the mesh. [`Self::garbage_collection`] must be
    /// called to remove the elements marked as deleted.
    pub fn delete_vertex(&mut self, delete_isolated_vertices: bool, v: VH) -> Result<(), Error> {
        self.topol
            .delete_vertex(v, delete_isolated_vertices, &mut self.cache)
    }

    /// Delete an edge.
    ///
    /// This will automatically delete all incident faces. Vertices in the
    /// neighborhood may become isolated. They are deleted if
    /// `delete_isolated_vertices` is true. The elements are marked as deleted,
    /// but not removed from the mesh. [`Self::garbage_collection`] must be
    /// called to remove the elements marked as deleted.
    pub fn delete_edge(&mut self, e: EH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_edge(
            e,
            delete_isolated_vertices,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    /// Delete a face.
    ///
    /// Vertices in the neighborhood may become isolated. They are deleted if
    /// `delete_isolated_vertices` is true. The elements are marked as deleted,
    /// but not removed from the mesh. [`Self::garbage_collection`] must be
    /// called to remove the elements marked as deleted.
    pub fn delete_face(&mut self, f: FH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_face(
            f,
            delete_isolated_vertices,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    /// Remove all elements in the mesh that are marked as deleted.
    ///
    /// The order of the elements may change during this operation. Because of
    /// this, all existing vertex handles may potentially become invalid, or
    /// point to the wrong elements. Properties defined on the mesh elements are
    /// safe because they are automatically synchronized.
    pub fn garbage_collection(&mut self) -> Result<(), Error> {
        self.topol.garbage_collection(&mut self.cache)
    }
}

impl<const DIM: usize, A> Clone for PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    /// Clone the mesh. This doesn't clone all the properties.
    ///
    /// This clones the topology of the mesh, i.e. the vertices, halfedges,
    /// edges, and faces, and all the built-in properties such as statuses of
    /// all the elements, points, and normals. This does not copy any other
    /// properties, because the mesh doesn't fully own the other properties. The
    /// properties are owned by whoever created them, and it is their
    /// responsibilities to properly clone the data in those
    /// properties. Furthermore, this only copies the data in built-in
    /// properties if they can be borrowed successfully. If someone upstream
    /// already borrowed these properties mutably, then the borrow during clone
    /// will fail and the data won't be copied. The new mesh will have default
    /// initialized properties.
    fn clone(&self) -> Self {
        let mut topol = self.topol.clone();
        let mut points = VProperty::new(&mut topol.vprops, A::zero_vector());
        match (self.points.try_borrow(), points.try_borrow_mut()) {
            (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
            _ => {}
        }
        let vnormals = match &self.vnormals {
            Some(normals) => {
                let mut dst = VProperty::new(&mut topol.vprops, A::zero_vector());
                match (normals.try_borrow(), dst.try_borrow_mut()) {
                    (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
                    _ => {}
                }
                Some(dst)
            }
            None => None,
        };
        let fnormals = match &self.fnormals {
            Some(normals) => {
                let mut dst = FProperty::new(&mut topol.fprops, A::zero_vector());
                match (normals.try_borrow(), dst.try_borrow_mut()) {
                    (Ok(src), Ok(mut dst)) if src.len() == dst.len() => dst.copy_from_slice(&src),
                    _ => {}
                }
                Some(dst)
            }
            None => None,
        };
        Self {
            topol,
            cache: Default::default(),
            points,
            vnormals,
            fnormals,
        }
    }
}

impl<const DIM: usize, A> HasTopology for PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    fn topology(&self) -> &Topology {
        &self.topol
    }

    fn topology_mut(&mut self) -> &mut Topology {
        &mut self.topol
    }
}

#[cfg(test)]
mod test {
    use crate::{use_glam::PolyMeshF32, Handle, HasTopology};

    #[test]
    fn t_icosahedron_clone() {
        let mut mesh = PolyMeshF32::icosahedron(1.0).expect("Cannot make an icosahedron");
        let myprop = mesh.create_halfedge_prop(0u8); // Create a custom property.
        {
            let myprop = myprop.try_borrow().expect("Cannot borrow property");
            assert_eq!(myprop.len(), mesh.num_halfedges());
        }
        // Tag the odd numbered elements.
        for v in mesh.vertices().filter(|v| v.index() % 2 != 0) {
            mesh.vertex_status_mut(v)
                .expect("Cannot access vertex status")
                .set_tagged(true);
        }
        for h in mesh.halfedges().filter(|h| h.index() % 2 != 0) {
            mesh.halfedge_status_mut(h)
                .expect("Cannot access halfedge status")
                .set_tagged(true);
        }
        for e in mesh.edges().filter(|e| e.index() % 2 != 0) {
            mesh.edge_status_mut(e)
                .expect("Cannot access edge status")
                .set_tagged(true);
        }
        for f in mesh.faces().filter(|f| f.index() % 2 != 0) {
            mesh.face_status_mut(f)
                .expect("Cannot access face status")
                .set_tagged(true);
        }
        let copy = mesh.clone();
        let src = mesh.points.try_borrow().expect("Cannot borrow points");
        let dst = copy.points.try_borrow().expect("Cannot borrow points");
        let src: &[glam::Vec3] = &src;
        let dst: &[glam::Vec3] = &dst;
        assert_eq!(src, dst);
        for v in copy.vertices() {
            assert_eq!(
                v.index() % 2 != 0,
                copy.vertex_status(v)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for h in copy.halfedges() {
            assert_eq!(
                h.index() % 2 != 0,
                copy.halfedge_status(h)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for e in copy.edges() {
            assert_eq!(
                e.index() % 2 != 0,
                copy.edge_status(e)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        for f in copy.faces() {
            assert_eq!(
                f.index() % 2 != 0,
                copy.face_status(f)
                    .expect("Cannot access vertex status")
                    .tagged()
            );
        }
        assert_eq!(
            mesh.try_calc_area().expect("Cannot compute area"),
            copy.try_calc_area().expect("Cannot compute area")
        );
        assert_eq!(
            mesh.try_calc_volume().expect("Cannot compute volume"),
            copy.try_calc_volume().expect("Cannot compute volume")
        );
        // Caller owned properties are not cloned. There is exactly one caller owned property.
        assert_eq!(1 + copy.num_halfedge_props(), mesh.num_halfedge_props());
    }
}
