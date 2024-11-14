use crate::{
    element::{EH, FH, HH, VH},
    error::Error,
    property::{EProperty, FProperty, HProperty, VProperty},
    status::Status,
    topol::{TopolCache, Topology},
};
use std::ops::Range;

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
        let points = topol.new_vprop(A::zero_vector());
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
        let points = topol.new_vprop_with_capacity(nverts, A::zero_vector());
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    /// Create a new vertex property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::alum_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_vertex_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_vertices(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    pub fn create_vertex_prop<T>(&mut self, default: T) -> VProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        self.topol.new_vprop(default)
    }

    /// Create a new halfedge property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::alum_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_halfedge_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_halfedges(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    pub fn create_halfedge_prop<T>(&mut self, default: T) -> HProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        self.topol.new_hprop(default)
    }

    /// Create a new edge property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::alum_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_edge_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_edges(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    pub fn create_edge_prop<T>(&mut self, default: T) -> EProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        self.topol.new_eprop(default)
    }

    /// Create a new face property of type T, with the `default` value.
    ///
    /// The default value will be used when new elements are added to the mesh.
    ///
    /// ```rust
    /// use alum::alum_glam::PolyMeshF32;
    ///
    /// let mut mesh = PolyMeshF32::tetrahedron(1.0).expect("Cannot crate tetrahedron");
    /// let prop = mesh.create_face_prop(42usize);
    /// // To use th property, you have to borrow it.
    /// let prop = prop.try_borrow().expect("Cannot borrow property");
    /// assert_eq!(mesh.num_faces(), prop.len());
    /// assert!(prop.iter().all(|v| *v == 42));
    /// ```
    pub fn create_face_prop<T>(&mut self, default: T) -> FProperty<T>
    where
        T: Clone + Copy + 'static,
    {
        self.topol.new_fprop(default)
    }

    /// Reserve memory for the given number of elements.
    ///
    /// The memory is also reserved for all properties.
    pub fn reserve(&mut self, nverts: usize, nedges: usize, nfaces: usize) -> Result<(), Error> {
        self.topol.reserve(nverts, nedges, nfaces)
    }

    /// Delete all elements and their properties.
    pub fn clear(&mut self) -> Result<(), Error> {
        self.topol.clear()
    }

    /// Number of vertices.
    pub fn num_vertices(&self) -> usize {
        self.topol.num_vertices()
    }

    /// Number of edges.
    pub fn num_edges(&self) -> usize {
        self.topol.num_edges()
    }

    /// Number of halfedges.
    pub fn num_halfedges(&self) -> usize {
        self.topol.num_halfedges()
    }

    /// Number of faces.
    pub fn num_faces(&self) -> usize {
        self.topol.num_faces()
    }

    /// Iterator over the vertices of the mesh.
    pub fn vertices(&self) -> impl Iterator<Item = VH> {
        self.topol.vertices()
    }

    /// Iterator over the halfedges of the mesh.
    pub fn halfedges(&self) -> impl Iterator<Item = HH> {
        self.topol.halfedges()
    }

    /// Iterator over the edges of the mesh.
    pub fn edges(&self) -> impl Iterator<Item = EH> {
        self.topol.edges()
    }

    /// Iterator over the faces of the mesh.
    pub fn faces(&self) -> impl Iterator<Item = FH> {
        self.topol.faces()
    }

    /// Check if this vertex is valid for this mesh.
    ///
    /// It's index has to be less than the number of vertices.
    pub fn is_valid_vertex(&self, v: VH) -> bool {
        self.topol.is_valid_vertex(v)
    }

    /// Check if this halfedge is valid for this mesh.
    ///
    /// It's index has to be less than the number of halfedges.
    pub fn is_valid_halfedge(&self, h: HH) -> bool {
        self.topol.is_valid_halfedge(h)
    }

    /// Check if this edge is valid for this mesh.
    ///
    /// It's index has to be less than the number of edges.
    pub fn is_valid_edge(&self, e: EH) -> bool {
        self.topol.is_valid_edge(e)
    }

    /// Check if the face is valid for this mesh.
    ///
    /// It's index has to be less than the number of faces.
    pub fn is_valid_face(&self, f: FH) -> bool {
        self.topol.is_valid_face(f)
    }

    /// Check if the vertex is manifold.
    ///
    /// A vertex is manifold if it has at most 1 outgoing halfedge.
    /// ```text
    ///    .......|     .......|.......     ....\     /...
    ///    .......|     .......|.......     .....\   /....
    ///    .......|     .......|.......     ......\ /.....
    ///    -------v     -------v-------     -------v------
    ///    .......|     .......|.......     ....../ \.....
    ///    .......|     .......|.......     ...../   \....
    ///    .......|     .......|.......     ..../     \...
    ///    Manifold     Manifold            Not manifold
    /// ```
    pub fn is_manifold_vertex(&self, v: VH) -> bool {
        self.topol.is_manifold_vertex(v)
    }

    /// Check if the vertex is on the boundary of the mesh.
    pub fn is_boundary_vertex(&self, v: VH) -> bool {
        self.topol.is_boundary_vertex(v)
    }

    /// Check if the halfedge is on the boundary.
    ///
    /// A halfedge is considered interior if it has a face incident on it.
    pub fn is_boundary_halfedge(&self, h: HH) -> bool {
        self.topol.is_boundary_halfedge(h)
    }

    /// Check if the edge is a boundary edge.
    ///
    /// An edge is considered interior if it has two faces incident on both of it's halfedges.
    pub fn is_boundary_edge(&self, e: EH) -> bool {
        self.topol.is_boundary_edge(e)
    }

    /// The number of edges incident on a vertex.
    pub fn vertex_valence(&self, v: VH) -> usize {
        self.topol.vertex_valence(v)
    }

    /// The number of vertices incident on a face.
    pub fn face_valence(&self, f: FH) -> usize {
        self.topol.face_valence(f)
    }

    /// The position of a vertex.
    pub fn point(&self, vi: VH) -> Result<A::Vector, Error> {
        self.points.get(vi)
    }

    /// The status of a vertex.
    pub fn vertex_status(&self, v: VH) -> Result<Status, Error> {
        self.topol.vertex_status(v)
    }

    /// The status of a halfedge.
    pub fn halfedge_status(&self, h: HH) -> Result<Status, Error> {
        self.topol.halfedge_status(h)
    }

    /// The status of an edge.
    pub fn edge_status(&self, e: EH) -> Result<Status, Error> {
        self.topol.edge_status(e)
    }

    /// The status of a face.
    pub fn face_status(&self, f: FH) -> Result<Status, Error> {
        self.topol.face_status(f)
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
    /// use alum::alum_glam::PolyMeshF32;
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
            .get_or_insert_with(|| self.topol.new_vprop(A::zero_vector()))
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
    /// use alum::alum_glam::PolyMeshF32;
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
            .get_or_insert_with(|| self.topol.new_fprop(A::zero_vector()))
            .clone()
    }

    /// Get the vertex the halfedge points to.
    pub fn head_vertex(&self, h: HH) -> VH {
        self.topol.head_vertex(h)
    }

    /// Get the vertex the halfedge is pointing away from.
    pub fn tail_vertex(&self, h: HH) -> VH {
        self.topol.tail_vertex(h)
    }

    /// Get the next halfedge in the loop.
    pub fn next_halfedge(&self, h: HH) -> HH {
        self.topol.next_halfedge(h)
    }

    /// Get the previous halfedge in the loop.
    pub fn prev_halfedge(&self, h: HH) -> HH {
        self.topol.prev_halfedge(h)
    }

    /// Get the opposite halfedge.
    pub fn opposite_halfedge(&self, h: HH) -> HH {
        self.topol.opposite_halfedge(h)
    }

    /// Get the face incident on the halfedge.
    pub fn halfedge_face(&self, h: HH) -> Option<FH> {
        self.topol.halfedge_face(h)
    }

    /// Get the edge corresponding to the halfedge.
    pub fn halfedge_edge(&self, h: HH) -> EH {
        self.topol.halfedge_edge(h)
    }

    /// Get the halfedge corresponding to the face.
    pub fn face_halfedge(&self, f: FH) -> HH {
        self.topol.face_halfedge(f)
    }

    /// Get the pair of halfedges associated with the given edge.
    pub fn halfedge_pair(&self, e: EH) -> (HH, HH) {
        self.topol.halfedge_pair(e)
    }

    /// Get a halfedge from the edge.
    ///
    /// The Boolean flag indicates one of the two possible orientations.
    pub fn edge_halfedge(&self, e: EH, flag: bool) -> HH {
        self.topol.edge_halfedge(e, flag)
    }

    /// Get the clockwise rotated halfedge around the vertex at the base of the
    /// given halfedge.
    pub fn cw_rotated_halfedge(&self, h: HH) -> HH {
        self.topol.cw_rotated_halfedge(h)
    }

    /// Get the counter-clockwise rotated hafedge around teh vertex at the base
    /// of the given halfedge.
    pub fn ccw_rotated_halfedge(&self, h: HH) -> HH {
        self.topol.ccw_rotated_halfedge(h)
    }

    /// Find a halfedge spanning the vertices `from` and `to`, if one exists
    pub fn find_halfedge(&self, from: VH, to: VH) -> Option<HH> {
        self.topol.find_halfedge(from, to)
    }

    /// Iterator over the vertex triplets that represent a triangulation of this
    /// mesh. The triangulation of a face does not take it's shape into
    /// account. It only accounts for the topology.
    ///
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, Handle};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_quad_face(0.into(), 1.into(), 2.into(), 3.into());
    /// assert_eq!(mesh.triangulated_vertices()
    ///                .flatten()
    ///                .map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 0, 1, 3, 1, 2]);
    /// ```
    pub fn triangulated_vertices(&self) -> impl Iterator<Item = [VH; 3]> + use<'_, A, DIM> {
        self.faces()
            .flat_map(move |f| self.triangulated_face_vertices(f))
    }

    /// Iterator over the vertex triplets that represent a triangulation of the
    /// given face.
    ///
    /// The triangulation does not take the shape of the face into account. It
    /// only accounts for the topology of the face.
    /// ```rust
    /// use alum::{alum_glam::PolyMeshF32, Handle};
    ///
    /// let mut mesh = PolyMeshF32::new();
    /// let verts = [glam::vec3(0.0, 0.0, 0.0), glam::vec3(1.0, 0.0, 0.0),
    ///              glam::vec3(1.0, 1.0, 0.0), glam::vec3(0.0, 1.0, 0.0)];
    /// mesh.add_vertices(&verts).expect("Cannot add vertices");
    /// mesh.add_quad_face(0.into(), 1.into(), 2.into(), 3.into());
    /// assert_eq!(mesh.triangulated_face_vertices(0.into())
    ///                .flatten()
    ///                .map(|v| v.index())
    ///                .collect::<Vec<u32>>(), [3, 0, 1, 3, 1, 2]);
    /// ```
    pub fn triangulated_face_vertices(
        &self,
        f: FH,
    ) -> impl Iterator<Item = [VH; 3]> + use<'_, A, DIM> {
        self.topol.triangulated_face_vertices(f)
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
    /// use alum::alum_glam::PolyMeshF32;
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
    pub fn add_vertices(&mut self, pos: &[A::Vector]) -> Result<Range<u32>, Error> {
        let vis = self.topol.add_vertices(pos.len())?;
        let mut points = self.points.try_borrow_mut()?;
        let points: &mut [A::Vector] = &mut points;
        for (i, v) in vis.clone().enumerate() {
            points[v as usize] = pos[i];
        }
        Ok(vis)
    }

    /// Add a face to this mesh. The `verts` are expected to be in counter-clockwise order.
    pub fn add_face(&mut self, verts: &[VH]) -> Result<FH, Error> {
        self.topol.add_face(verts, &mut self.cache)
    }

    /// Add a triangular face with the three given vertices.
    pub fn add_tri_face(&mut self, v0: VH, v1: VH, v2: VH) -> Result<FH, Error> {
        self.add_face(&[v0, v1, v2])
    }

    /// Add a quadrilateral face with the three given vertices.
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
