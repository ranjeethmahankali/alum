pub mod error;
pub mod math;
pub mod mesh;
pub mod obj;
pub mod property;
pub mod status;
pub mod topol;

mod edit;
mod element;
mod iterator;
mod macros;
mod primitive;

#[cfg(feature = "use_glam")]
mod use_glam;

#[cfg(feature = "use_glam")]
pub type PolyMeshF32 = mesh::PolyMeshT<3, use_glam::BuiltInAdaptorF32>;

#[cfg(feature = "use_glam")]
pub type PolyMeshF64 = mesh::PolyMeshT<3, use_glam::BuiltInAdaptorF64>;

pub use element::{Handle, EH, FH, HH, VH};
pub use mesh::PolyMeshT;
pub use topol::Topology;
