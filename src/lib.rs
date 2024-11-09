pub mod element;
pub mod error;
pub mod mesh;
pub mod obj;
pub mod primitive;
pub mod property;
pub mod status;
pub mod topol;
pub mod vector;

mod iterator;
mod macros;
mod triangulate;
mod collapse;

pub type PolyMeshF32 = mesh::PolyMeshF32;
pub type Topology = topol::Topology;
