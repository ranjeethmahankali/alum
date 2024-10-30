use crate::element::{FH, HH, VH};

#[derive(Debug)]
pub enum Error {
    // Properties.
    BorrowedPropertyAccess,
    PropertyDoesNotExist,
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
}
