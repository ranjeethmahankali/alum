use std::fmt::{Debug, Display};

pub trait Handle {
    fn index(&self) -> u32;
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct VH {
    idx: u32,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct HH {
    idx: u32,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct EH {
    idx: u32,
}

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