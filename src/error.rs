use std::path::PathBuf;

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
    // Obj.
    InvalidObjFile(PathBuf),
    ObjLoadFailed(String),
    IncorrectNumberOfCoordinates(usize),
    IncorrectIndexCount(usize),
    // Other,
    MismatchedArrayLengths(usize, usize),
}
