use std::path::PathBuf;

use crate::element::{FH, HH, VH};

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
    /// Trying to access a face that has been deleted.
    DeletedFace(FH),
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
    /// The mesh has deleted elements, which create problems for further
    /// operations. Garbage collection needs to be performed before proceeding.
    GarbageCollectionRequired,
}
