use std::{fmt::Debug, path::PathBuf};

use crate::{
    element::{FH, HH, VH},
    EH,
};

/// Error type for all operations across this crate.
#[derive(Debug)]
pub enum Error {
    /// Trying to borrow a property that is already borrowed.
    BorrowedPropertyAccess,
    /// Accessing a property that has been dropped.
    PropertyDoesNotExist,
    /// Face have not been computed for the given mesh.
    FaceNormalsNotAvailable,
    /// Trying to access outside the bounds of a property vector.
    OutOfBoundsAccess,
    /// Non manifold topology at a vertex.
    ComplexVertex(VH),
    /// Non manifold topology at a halfedge.
    ComplexHalfedge(HH),
    /// Failed to create valid topology when adding a face.
    PatchRelinkingFailed,
    /// Encountered a face that has been deleted.
    DeletedFace(FH),
    /// Encountered an edge that has been deleted.
    DeletedEdge(EH),
    /// Encountered a halfedge that has been deleted.
    DeletedHalfedge(HH),
    /// Encountered a vertex that has been deleted.
    DeletedVertex(VH),
    /// The given OBJ file is not valid.
    InvalidObjFile(PathBuf),
    /// Failed to load a mesh from the given obj file.
    ObjLoadFailed(String),
    /// The given is of the wrong length to represent the coordinates of 3d
    /// points.
    IncorrectNumberOfCoordinates(usize),
    /// Two arrays whose lengths are expected to be the same, have different
    /// lengths.
    MismatchedArrayLengths(usize, usize),

    // Topological checks.
    /// The outgoing halfedge of a boundary vertex must be a boundary halfedge.
    OutgoingHalfedgeNotBoundary(VH),
    /// Invalid halfedge.
    InvalidHalfedge(HH),
    /// Broken topology found among the outgoing halfedges around a vertex.
    InvalidOutgoingHalfedges(VH),
    /// Halfedge has the same tail and head vertex.
    DegenerateHalfedge(HH),

    /// Cycle found when iterating over halfedges in a face loop.
    InvalidLoopTopology(HH),
    /// Halfedges in a loop are not linked to the same face / boundary.
    InconsistentFaceInLoop(HH),
    /// The next-previous link between a halfedge pair doesn't commute.
    InvalidHalfedgeLink(HH),
    /// The link between halfedge and vertex is inconsistent.
    InvalidHalfedgeVertexLink(HH),
    /// Invalid face-halfedge link.
    InvalidFaceHalfedgeLink(FH, HH),

    /// Edge is not a unique link between its incident faces. Removing it will
    /// produce dangling / invalid topology.
    EdgeIsNotAUniqueLink(EH),
    /// Removing a boundary edge breaks the topology of the mesh. If you really
    /// want to do this, considering deleting the edge.
    CannotRemoveBoundaryEdge(EH),
    /// The halfedges are expected to be in the same loop, but they are not.
    HalfedgesNotInTheSameLoop(HH, HH),
    /// Cannot insert an edge spanning the two halfedges.
    CannotInsertEdge(HH, HH),

    /// Something went wrong when trying to split a face.
    CannotSplitFace(FH),

    /// Cannot swap a boundary edge.
    CannotSwapBoundaryEdge(EH),
    /// Degenerate edge.
    DegenerateEdge(EH),
}
