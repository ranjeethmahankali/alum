use crate::topol::{HH, VH};

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
}
