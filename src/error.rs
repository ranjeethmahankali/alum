use crate::element::{FH, HH, VH};

#[derive(Debug)]
pub enum Error {
    // Properties.
    BorrowedPropertyAccess,
    PropertyDoesNotExist,
    FaceNormalsNotAvailable,
    OutOfBoundsAccess,
    // Topology.
    ComplexVertex(VH),
    ComplexHalfedge(HH),
    HalfedgeNotFound,
    PatchRelinkingFailed,
    DeletedFace(FH),
    MismatchedArrayLengths(usize, usize),
    // Obj.
    ObjLoadFailed(String),
    IncorrectNumberOfCoordinates(usize),
    IncorrectIndexCount(usize),
}
