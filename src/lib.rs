pub mod element;
pub mod error;
pub mod mesh;
pub mod obj;
pub mod primitive;
pub mod property;
pub mod status;
pub mod topol;
pub mod triangulate;
pub mod vector;

mod collapse;
mod iterator;
mod macros;

pub type PolyMeshF32 = mesh::PolyMeshF32;
pub type Topology = topol::Topology;
