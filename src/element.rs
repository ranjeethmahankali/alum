use crate::{iterator, topol::Topology, Adaptor, PolyMeshT};
use std::fmt::{Debug, Display};

/**
 * All elements of the mesh implement this trait. They are identified by their
 * index.
 */
pub trait Handle {
    /**
     * The index of the element.
     */
    fn index(&self) -> u32;
}

/**
 * Vertex handle.
 */
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VH {
    idx: u32,
}

/**
 * Halfedge handle.
 */
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct HH {
    idx: u32,
}

/**
 * Edge handle.
 */
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct EH {
    idx: u32,
}

/**
 * Face handle.
 */
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

pub trait HasTopology {
    fn topology(&self) -> &Topology;
}

impl VH {
    pub fn halfedge(self, mesh: &impl HasTopology) -> Option<HH> {
        mesh.topology().vertex_halfedge(self)
    }

    /// Check if this vertex is valid for the `mesh`.
    ///
    /// The index has to be less than the number of vertices in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_valid_vertex(self)
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
    pub fn is_manifold(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_manifold_vertex(self)
    }

    /// Check if this vertex is on the boundary of the `mesh`.
    pub fn is_boundary(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_boundary_vertex(self)
    }

    /// The number of edges incident on this vertex.
    pub fn valence(self, mesh: &impl HasTopology) -> usize {
        iterator::voh_ccw_iter(mesh.topology(), self).count()
    }
}

impl HH {
    pub fn head(self, mesh: &impl HasTopology) -> VH {
        mesh.topology().head_vertex(self)
    }

    pub fn tail(self, mesh: &impl HasTopology) -> VH {
        mesh.topology().tail_vertex(self)
    }

    pub fn opposite(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().opposite_halfedge(self)
    }

    pub fn prev(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().prev_halfedge(self)
    }

    pub fn next(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().next_halfedge(self)
    }

    pub fn face(self, mesh: &impl HasTopology) -> Option<FH> {
        mesh.topology().halfedge_face(self)
    }

    /// Check if this halfedge is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halfedges in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_valid_halfedge(self)
    }

    /// Check if this halfedge is on the boundary of `mesh`.
    ///
    /// A halfedge is considered interior if it has a face incident on it.
    pub fn is_boundary(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_boundary_halfedge(self)
    }
}

impl EH {
    pub fn halfedges(self) -> (HH, HH) {
        let hi = self.idx << 1;
        (hi.into(), (hi | 1).into())
    }

    pub fn halfedge(self, flag: bool) -> HH {
        ((self.idx << 1) | if flag { 1 } else { 0 }).into()
    }

    /// Check if this edge is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halfedges in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_valid_edge(self)
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
    pub fn halfedge(self, mesh: &impl HasTopology) -> HH {
        mesh.topology().face_halfedge(self)
    }

    /// Check if this face is valid for the `mesh`.
    ///
    /// The index has to be less than the number of halffaces in the mesh.
    pub fn is_valid(self, mesh: &impl HasTopology) -> bool {
        mesh.topology().is_valid_face(self)
    }
}

impl HasTopology for Topology {
    fn topology(&self) -> &Topology {
        self
    }
}

impl<const DIM: usize, A> HasTopology for PolyMeshT<DIM, A>
where
    A: Adaptor<DIM>,
{
    fn topology(&self) -> &Topology {
        &self.topol
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
