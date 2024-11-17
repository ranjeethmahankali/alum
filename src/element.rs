use crate::{iterator::HasIterators, topol::HasTopology};
use std::{
    fmt::{Debug, Display},
    marker::PhantomData,
    ops::Range,
};

/// All elements of the mesh implement this trait. They are identified by their
/// index.
pub trait Handle: From<u32> {
    /// The index of the element.
    fn index(&self) -> u32;
}

pub struct HandleRange<H>
where
    H: Handle,
{
    current: u32,
    stop: u32,
    _phantom: PhantomData<H>,
}

impl<H> From<Range<u32>> for HandleRange<H>
where
    H: Handle,
{
    fn from(value: Range<u32>) -> Self {
        HandleRange {
            current: value.start,
            stop: value.end,
            _phantom: PhantomData,
        }
    }
}

impl<H> Iterator for HandleRange<H>
where
    H: Handle,
{
    type Item = H;

    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.stop {
            let i = self.current;
            self.current += 1;
            Some(i.into())
        } else {
            None
        }
    }
}

pub type VRange = HandleRange<VH>;
pub type HRange = HandleRange<HH>;
pub type ERange = HandleRange<EH>;
pub type FRange = HandleRange<FH>;

/// Vertex handle.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VH {
    idx: u32,
}

/// Halfedge handle.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct HH {
    idx: u32,
}

/// Edge handle.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct EH {
    idx: u32,
}

/// Face handle.
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct FH {
    idx: u32,
}

impl Handle for VH {
    fn index(&self) -> u32 {
        self.idx
    }
}

impl From<u32> for VH {
    fn from(idx: u32) -> Self {
        VH { idx }
    }
}

impl From<&u32> for VH {
    fn from(idx: &u32) -> Self {
        VH { idx: *idx }
    }
}

impl From<&mut u32> for VH {
    fn from(idx: &mut u32) -> Self {
        VH { idx: *idx }
    }
}

impl Handle for HH {
    fn index(&self) -> u32 {
        self.idx
    }
}

impl From<u32> for HH {
    fn from(idx: u32) -> Self {
        HH { idx }
    }
}

impl From<&u32> for HH {
    fn from(idx: &u32) -> Self {
        HH { idx: *idx }
    }
}

impl From<&mut u32> for HH {
    fn from(idx: &mut u32) -> Self {
        HH { idx: *idx }
    }
}

impl Handle for EH {
    fn index(&self) -> u32 {
        self.idx
    }
}

impl From<u32> for EH {
    fn from(idx: u32) -> Self {
        EH { idx }
    }
}

impl From<&mut u32> for EH {
    fn from(idx: &mut u32) -> Self {
        EH { idx: *idx }
    }
}

impl Handle for FH {
    fn index(&self) -> u32 {
        self.idx
    }
}

impl From<u32> for FH {
    fn from(idx: u32) -> Self {
        FH { idx }
    }
}

impl From<&u32> for FH {
    fn from(idx: &u32) -> Self {
        FH { idx: *idx }
    }
}

impl From<&mut u32> for FH {
    fn from(idx: &mut u32) -> Self {
        FH { idx: *idx }
    }
}

impl Display for VH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "VH({})", self.index())
    }
}

impl Display for HH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HH({})", self.index())
    }
}

impl Display for EH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "EH({})", self.index())
    }
}

impl Display for FH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FH({})", self.index())
    }
}

impl Debug for VH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "VH({})", self.index())
    }
}

impl Debug for HH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "HH({})", self.index())
    }
}

impl Debug for EH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "EH({})", self.index())
    }
}

impl Debug for FH {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "FH({})", self.index())
    }
}

impl VH {
    pub fn halfedge(self, mesh: &impl HasTopology) -> Option<HH> {
        mesh.topology().vertex(self).halfedge
    }

    /// Check if this vertex is valid for the `mesh`.
    ///
    /// The index has to be less than the number of vertices in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        (self.idx as usize) < mesh.topology().num_vertices()
    }

    /// Check if this vertex is manifold.
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
    pub fn is_manifold(self, mesh: &impl HasIterators) -> bool {
        mesh.voh_ccw_iter(self)
            .skip(1)
            .all(|h| !h.is_boundary(mesh))
    }

    /// Check if this vertex is on the boundary of the `mesh`.
    pub fn is_boundary(self, mesh: &impl HasTopology) -> bool {
        match self.halfedge(mesh) {
            Some(h) => h.is_boundary(mesh),
            None => true,
        }
    }

    /// The number of edges incident on this vertex.
    pub fn valence(self, mesh: &impl HasIterators) -> usize {
        mesh.voh_ccw_iter(self).count()
    }
}

impl HH {
    /// Get the vertex the halfedge points to.
    pub fn head(self, mesh: &impl HasTopology) -> VH {
        mesh.topology().halfedge(self).vertex
    }

    /// Get the vertex the halfedge is pointing away from.
    pub fn tail(self, mesh: &impl HasTopology) -> VH {
        mesh.topology().halfedge(self.opposite()).vertex
    }

    /// Get the opposite halfedge.
    pub fn opposite(self) -> HH {
        (self.idx ^ 1).into()
    }

    /// Get the previous halfedge in the loop.
    pub fn prev(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().halfedge(self).prev
    }

    /// Get the next halfedge in the loop.
    pub fn next(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().halfedge(self).next
    }

    /// Get the face incident on the halfedge.
    pub fn face(self, mesh: &impl HasTopology) -> Option<FH> {
        mesh.topology().halfedge(self).face
    }

    /// Get the edge corresponding to the halfedge.
    pub fn edge(self) -> EH {
        (self.idx >> 1).into()
    }

    /// Check if this halfedge is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halfedges in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        (self.idx as usize) < mesh.topology().num_halfedges()
    }

    /// Check if this halfedge is on the boundary of `mesh`.
    ///
    /// A halfedge is considered interior if it has a face incident on it.
    pub fn is_boundary(self, mesh: &impl HasTopology) -> bool {
        self.face(mesh).is_none()
    }

    /// Get the clockwise rotated halfedge around the vertex at the base of the
    /// given halfedge.
    pub fn rotate_cw(self, mesh: &impl HasTopology) -> HH {
        self.opposite().next(mesh)
    }

    /// Get the counter-clockwise rotated hafedge around the vertex at the base
    /// of the given halfedge.
    pub fn rotate_ccw(self, mesh: &impl HasTopology) -> HH {
        self.prev(mesh).opposite()
    }
}

impl EH {
    /// Get the pair of halfedges associated with the given edge.
    pub fn halfedges(self) -> (HH, HH) {
        let hi = self.idx << 1;
        (hi.into(), (hi | 1).into())
    }

    /// Get the two vertices incident on this edge.
    pub fn vertices(self, mesh: &impl HasTopology) -> (VH, VH) {
        let (h, oh) = self.halfedges();
        (h.head(mesh), oh.head(mesh))
    }

    /// Get a halfedge from the edge.
    ///
    /// The Boolean flag indicates one of the two possible orientations.
    pub fn halfedge(self, flag: bool) -> HH {
        ((self.idx << 1) | if flag { 1 } else { 0 }).into()
    }

    /// Check if this edge is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halfedges in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        (self.idx as usize) < mesh.topology().num_edges()
    }

    /// Check if the edge is a boundary edge.
    ///
    /// An edge is considered interior if it has two faces incident on both of it's halfedges.
    pub fn is_boundary(self, mesh: &impl HasTopology) -> bool {
        let (h, oh) = self.halfedges();
        h.is_boundary(mesh) || oh.is_boundary(mesh)
    }
}

impl FH {
    /// Get the halfedge corresponding to the face.
    pub fn halfedge(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().face(self).halfedge
    }

    /// Check if this face is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halffaces in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        (self.idx as usize) < mesh.topology().num_faces()
    }

    /// The number of vertices incident on a face.
    pub fn valence(self, mesh: &impl HasIterators) -> usize {
        mesh.fh_ccw_iter(self).count()
    }

    /// Check if a face is on the boundary.
    ///
    /// A face is considered to be on the boundary if at least one of it's edges
    /// are on the boundary. If `check_vertex` is true, then it is also
    /// considered to be on the boundary if at least one incident vertex is on
    /// the boundary.
    pub fn is_boundary(self, mesh: &impl HasIterators, check_vertex: bool) -> bool {
        mesh.fh_ccw_iter(self)
            .any(|h| h.opposite().is_boundary(mesh))
            || (check_vertex && mesh.fv_ccw_iter(self).any(|v| v.is_boundary(mesh)))
    }
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct Vertex {
    pub(crate) halfedge: Option<HH>,
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct Halfedge {
    pub(crate) face: Option<FH>,
    pub(crate) vertex: VH,
    pub(crate) next: HH,
    pub(crate) prev: HH,
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct Edge {
    pub(crate) halfedges: [Halfedge; 2],
}

#[derive(Debug, Copy, Clone)]
pub(crate) struct Face {
    pub(crate) halfedge: HH,
}
