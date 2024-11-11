use crate::{
    element::{Handle, EH, FH, HH, VH},
    error::Error,
    iterator,
    property::{EProperty, FProperty, HProperty, VProperty},
    status::Status,
    topol::{TopolCache, Topology},
    vector::TVec,
};

/// A polygon mesh. `VecT` is the 3d vector type used to represent this
/// mesh. The positions of the vertices and face normals etc. are all
/// `VecT`. Quantities like edge length, face area etc. are of type
/// `VecT::Scalar`. Various functions of , such as computing normals
/// etc. are available predicated on the trait bounds `VecT` satisfies.
pub struct PolyMeshT<VecT, const DIM: usize>
where
    VecT: TVec<DIM>,
{
    pub(crate) topol: Topology,
    pub(crate) cache: TopolCache,
    points: VProperty<VecT>,
    vnormals: Option<VProperty<VecT>>,
    fnormals: Option<FProperty<VecT>>,
}

pub type PolyMeshF32 = PolyMeshT<glam::Vec3, 3>;
pub type PolyMeshF64 = PolyMeshT<glam::DVec3, 3>;

impl<VecT, const DIM: usize> Default for PolyMeshT<VecT, DIM>
where
    VecT: TVec<DIM>,
{
    /// Create a new empty mesh.
    fn default() -> Self {
        Self::new()
    }
}

impl<VecT, const DIM: usize> PolyMeshT<VecT, DIM>
where
    VecT: TVec<DIM>,
{
    /// Create a new empty mesh.
    pub fn new() -> Self {
        let mut topol = Topology::new();
        let points = topol.new_vprop();
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    /// Create an empty mesh with memory reserved for the given number of
    /// elements.  The memory is also reserved for all the built-in
    /// properties. Normals are not computed.
    pub fn with_capacity(nverts: usize, nedges: usize, nfaces: usize) -> Self {
        let mut topol = Topology::with_capacity(nverts, nedges, nfaces);
        let points = topol.new_vprop_with_capacity(nverts);
        PolyMeshT {
            topol,
            cache: TopolCache::default(),
            points,
            vnormals: None,
            fnormals: None,
        }
    }

    /// Create a new vertex property of type T.
    pub fn create_vertex_prop<T>(&mut self) -> VProperty<T>
    where
        T: Default + Clone + Copy + 'static,
    {
        self.topol.new_vprop()
    }

    /// Create a new halfedge property of type T.
    pub fn create_halfedge_prop<T>(&mut self) -> HProperty<T>
    where
        T: Default + Clone + Copy + 'static,
    {
        self.topol.new_hprop()
    }

    /// Create a new edge property of type T.
    pub fn create_edge_prop<T>(&mut self) -> EProperty<T>
    where
        T: Default + Clone + Copy + 'static,
    {
        self.topol.new_eprop()
    }

    /// Create a new face property of type T.
    pub fn create_face_prop<T>(&mut self) -> FProperty<T>
    where
        T: Default + Clone + Copy + 'static,
    {
        self.topol.new_fprop()
    }

    /// Reserve memory for the given number of elements. The memory is also
    /// reserved for all properties.
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

    /// Get the topology of this mesh.
    pub fn topology(&self) -> &Topology {
        &self.topol
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

    /// Check if this vertex is valid for this mesh. It's index has to be less
    /// than the number of vertices.
    pub fn is_valid_vertex(&self, v: VH) -> bool {
        self.topol.is_valid_vertex(v)
    }

    /// Check if this halfedge is valid for this mesh. It's index has to be less
    /// than the number of halfedges.
    pub fn is_valid_halfedge(&self, h: HH) -> bool {
        self.topol.is_valid_halfedge(h)
    }

    /// Check if this edge is valid for this mesh. It's index has to be less
    /// than the number of edges.
    pub fn is_valid_edge(&self, e: EH) -> bool {
        self.topol.is_valid_edge(e)
    }

    /// Check if the face is valid for this mesh. It's index has to be less than
    /// the number of faces.
    pub fn is_valid_face(&self, f: FH) -> bool {
        self.topol.is_valid_face(f)
    }

    /// Check if the vertex is manifold. A vertex is manifold if it has at most 1 outgoing halfedge.
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

    /// Check if the halfedge is on the boundary. A halfedge is considered
    /// interior if it has a face incident on it.
    pub fn is_boundary_halfedge(&self, h: HH) -> bool {
        self.topol.is_boundary_halfedge(h)
    }

    /// Check if the edge is a boundary edge.
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
    pub fn point(&self, vi: VH) -> Result<VecT, Error> {
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
    pub fn set_point(&mut self, v: VH, pos: VecT) -> Result<(), Error> {
        self.points.set(v, pos)
    }

    /// Get the property corresponding to the positions of the vertices.
    pub fn points(&self) -> VProperty<VecT> {
        self.points.clone()
    }

    /// Check if the vertex normals are available.
    pub fn has_vertex_normals(&self) -> bool {
        self.vnormals.is_some()
    }

    /// Get the vertex normals. It is `None` if the vertex normals are not
    /// available.
    pub fn vertex_normals(&self) -> Option<VProperty<VecT>> {
        self.vnormals.clone()
    }

    /// If the vertex normals property is available, it is returned as is. If it
    /// is not available, a new property for vertex normals is created and
    /// returned. The vertex normals are not computed. The values in the
    /// property are default initialized.
    pub fn request_vertex_normals(&mut self) -> VProperty<VecT> {
        self.vnormals
            .get_or_insert_with(|| self.topol.new_vprop())
            .clone()
    }

    /// Check if the face normals are available.
    pub fn has_face_normals(&self) -> bool {
        self.fnormals.is_some()
    }

    /// Get the face normals. It is `None` if the face normals are not available.
    pub fn face_normals(&self) -> Option<FProperty<VecT>> {
        self.fnormals.clone()
    }

    /// If the face normals property is available, it is returned as is. If it is
    /// not available, a new property for face normals is created and
    /// returned. The face normals are not computed. The values in the property
    /// are default initialized.
    pub fn request_face_normals(&mut self) -> FProperty<VecT> {
        self.fnormals
            .get_or_insert_with(|| self.topol.new_fprop())
            .clone()
    }

    /// Get the vertex the halfedge points to.
    pub fn to_vertex(&self, h: HH) -> VH {
        self.topol.to_vertex(h)
    }

    /// Get the vertex the halfedge is pointing away from.
    pub fn from_vertex(&self, h: HH) -> VH {
        self.topol.from_vertex(h)
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

    /// Get a halfedge from the edge. The Boolean flag indicates one of the two
    /// possible orientations.
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

    /// Iterator over the outgoing halfedges around a vertex, going counter-clockwise.
    pub fn voh_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::voh_ccw_iter(&self.topol, v)
    }

    /// Iterator over the outgoing halfedges around a vertex, going clockwise
    pub fn voh_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::voh_cw_iter(&self.topol, v)
    }

    /// Iterator over the incoming halfedges around a vertex, going
    /// counter-clockwise
    pub fn vih_ccw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::vih_ccw_iter(&self.topol, v)
    }

    /// Iterator over the incoming halfedges around a vertex, going clockwise
    pub fn vih_cw_iter(&self, v: VH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::vih_cw_iter(&self.topol, v)
    }

    /// Iterator over the faces incident on a vertex, going counter-clockwise.
    pub fn vf_ccw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, VecT, DIM> {
        iterator::vf_ccw_iter(&self.topol, v)
    }

    /// Iterator over the faces incident on a vertex, going clockwise.
    pub fn vf_cw_iter(&self, v: VH) -> impl Iterator<Item = FH> + use<'_, VecT, DIM> {
        iterator::vf_cw_iter(&self.topol, v)
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// counter-clockwise.
    pub fn vv_ccw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, VecT, DIM> {
        iterator::vv_ccw_iter(&self.topol, v)
    }

    /// Iterator over the neighboring vertices around the given vertex, going
    /// clockwise.
    pub fn vv_cw_iter(&self, v: VH) -> impl Iterator<Item = VH> + use<'_, VecT, DIM> {
        iterator::vv_cw_iter(&self.topol, v)
    }

    /// Iterator over the incident edges around an vertex, going counter-clockwise.
    pub fn ve_ccw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, VecT, DIM> {
        iterator::ve_ccw_iter(&self.topol, v)
    }

    /// Iterator over the incident edges around a vertex, going clockwise.
    pub fn ve_cw_iter(&self, v: VH) -> impl Iterator<Item = EH> + use<'_, VecT, DIM> {
        iterator::ve_cw_iter(&self.topol, v)
    }

    /// Iterator over the two vertices incident on the given edge.
    pub fn ev_iter(&self, e: EH) -> impl Iterator<Item = VH> + use<'_, VecT, DIM> {
        iterator::ev_iter(&self.topol, e)
    }

    /// Iterator over the two halfedges corresponding to an edge.
    pub fn eh_iter(&self, e: EH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::eh_iter(&self.topol, e)
    }

    /// Iterator over the faces incident on an edge.
    pub fn ef_iter(&self, e: EH) -> impl Iterator<Item = FH> + use<'_, VecT, DIM> {
        iterator::ef_iter(&self.topol, e)
    }

    /// Iterator over the halfedges of a face loop, going counter-clockwise.
    pub fn fh_ccw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::fh_ccw_iter(&self.topol, f)
    }

    /// Iterator over the halfedges of a face loop, going clockwise.
    pub fn fh_cw_iter(&self, f: FH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::fh_cw_iter(&self.topol, f)
    }

    /// Iterator over the vertices incident on a face, going counter-clockwise.
    pub fn fv_ccw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, VecT, DIM> {
        iterator::fv_ccw_iter(&self.topol, f)
    }

    /// Iterator over the vertices incident on a face, going clockwise.
    pub fn fv_cw_iter(&self, f: FH) -> impl Iterator<Item = VH> + use<'_, VecT, DIM> {
        iterator::fv_cw_iter(&self.topol, f)
    }

    /// Iterator over the edges incident on a face, going counter-clockwise.
    pub fn fe_ccw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, VecT, DIM> {
        iterator::fe_ccw_iter(&self.topol, f)
    }

    /// Iterator over the edges incident on a face, going clockwise.
    pub fn fe_cw_iter(&self, f: FH) -> impl Iterator<Item = EH> + use<'_, VecT, DIM> {
        iterator::fe_cw_iter(&self.topol, f)
    }

    /// Iterator over the neighboring faces of the given face, going
    /// counter-clockwise. This includes the faces connected via a shared, edge,
    /// but not those connected via a shared vertex.
    pub fn ff_ccw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, VecT, DIM> {
        iterator::ff_ccw_iter(&self.topol, f)
    }

    /// Iterator over the neighboring faces of the given face, going
    /// clockwise. This includes the faces connected via a shared, edge, but not
    /// those connected via a shared vertex.
    pub fn ff_cw_iter(&self, f: FH) -> impl Iterator<Item = FH> + use<'_, VecT, DIM> {
        iterator::ff_cw_iter(&self.topol, f)
    }

    /// This is similar to `voh_ccw_iter` around the base of the given halfedge,
    /// except this iterator starts at the provided halfedge. So this is
    /// equivalent to a circular shifted `voh_ccw_iter` of the vertex at teh base
    /// of this halfedge.
    pub fn ccw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::ccw_rotate_iter(&self.topol, h)
    }

    /// This is similar to `voh_cw_iter` around the base of the given halfedge,
    /// except this iterator starts at the provided halfedge. So this is
    /// equivalent to a circular shifted `voh_cw_iter` of the vertex at teh base
    /// of this halfedge.
    pub fn cw_rotate_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::cw_rotate_iter(&self.topol, h)
    }

    /// Counter-clockwise iterator over the halfedges in a loop. The iterator
    /// will start at the given halfedge. If the halfedge has an incident face,
    /// this iterator is equivalent to a circular shifted `fh_ccw_iter` of the
    /// incident face. If the halfedge is on the boundary, this iterator goes
    /// over the boundary loop counter-clockwise.
    pub fn loop_ccw_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::loop_ccw_iter(&self.topol, h)
    }

    /// Counter-clockwise iterator over the halfedges in a loop. The iterator
    /// will start at the given halfedge. If the halfedge has an incident face,
    /// this iterator is equivalent to a circular shifted `fh_cw_iter` of the
    /// incident face. If the halfedge is on the boundary, this iterator goes
    /// over the boundary loop clockwise.
    pub fn loop_cw_iter(&self, h: HH) -> impl Iterator<Item = HH> + use<'_, VecT, DIM> {
        iterator::loop_cw_iter(&self.topol, h)
    }

    /// Iterator over the vertex triplets that represent a triangulation of this
    /// mesh. The triangulation of a face does not take it's shape into
    /// account. It only accounts for the topology.
    pub fn triangulated_vertices(&self) -> impl Iterator<Item = [VH; 3]> + use<'_, VecT, DIM> {
        self.topol.triangulated_vertices()
    }

    /// Iterator over the vertex triplets that represent a triangulation of the
    /// given face. The triangulation does not take the shape of the face into
    /// account. It only accounts for the topology of the face.
    pub fn triangulated_face_vertices(
        &self,
        f: FH,
    ) -> impl Iterator<Item = [VH; 3]> + use<'_, VecT, DIM> {
        self.topol.triangulated_face_vertices(f)
    }

    /// Add a vertex to this mesh at the given position.
    pub fn add_vertex(&mut self, pos: VecT) -> Result<VH, Error> {
        let vi = self.topol.add_vertex()?;
        self.points.set(vi, pos)?;
        Ok(vi)
    }

    /// Add several vertices at once. The positions of the vertices must be
    /// supplied in `pos`. The handles of the added vertices will be written into
    /// `dst`.
    pub fn add_vertices(&mut self, pos: &[VecT], dst: &mut [VH]) -> Result<(), Error> {
        if pos.len() != dst.len() {
            return Err(Error::MismatchedArrayLengths(pos.len(), dst.len()));
        }
        self.topol.add_vertices(dst)?;
        let verts: &[VH] = dst; // Make immutable.
        {
            let mut points = self.points.try_borrow_mut()?;
            let points: &mut [VecT] = &mut points;
            for (i, v) in verts.iter().enumerate() {
                points[v.index() as usize] = pos[i];
            }
        }
        Ok(())
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

    /// Delete a vertex. This will automatically delete all incident edges and
    /// faces. Vertices in the neighborhood may become isolated. They are deleted
    /// if `delete_isolated_vertices` is true. The elements are marked as
    /// deleted, but not removed from the mesh. `garbage_collection` must be
    /// called to remove the elements marked as deleted.
    pub fn delete_vertex(&mut self, delete_isolated_vertices: bool, v: VH) -> Result<(), Error> {
        self.topol
            .delete_vertex(v, delete_isolated_vertices, &mut self.cache)
    }

    /// Delete an edge. This will automatically delete all incident
    /// faces. Vertices in the neighborhood may become isolated. They are deleted
    /// if `delete_isolated_vertices` is true. The elements are marked as
    /// deleted, but not removed from the mesh. `garbage_collection` must be
    /// called to remove the elements marked as deleted.
    pub fn delete_edge(&mut self, e: EH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_edge(
            e,
            delete_isolated_vertices,
            &mut self.cache.halfedges,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    /// Delete a face. Vertices in the neighborhood may become isolated. They are
    /// deleted if `delete_isolated_vertices` is true. The elements are marked as
    /// deleted, but not removed from the mesh. `garbage_collection` must be
    /// called to remove the elements marked as deleted.
    pub fn delete_face(&mut self, f: FH, delete_isolated_vertices: bool) -> Result<(), Error> {
        self.topol.delete_face(
            f,
            delete_isolated_vertices,
            &mut self.cache.halfedges,
            &mut self.cache.edges,
            &mut self.cache.vertices,
        )
    }

    /// Remove all elements in the mesh that are marked as deleted. The order of
    /// the elements may change during this operation, but all the porperty
    /// vectors are synchronized.
    pub fn garbage_collection(&mut self) -> Result<(), Error> {
        self.topol.garbage_collection(&mut self.cache)
    }
}
