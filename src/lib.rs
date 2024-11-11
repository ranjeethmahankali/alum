pub mod adaptor;
pub mod element;
pub mod error;
pub mod math;
pub mod mesh;
pub mod obj;
pub mod primitive;
pub mod property;
pub mod status;
pub mod topol;

mod edit;
mod iterator;
mod macros;

pub type PolyMeshF32 = mesh::PolyMeshT<3, adaptor::DefaultAdaptorF32>;
pub type PolyMeshF64 = mesh::PolyMeshT<3, adaptor::DefaultAdaptorF64>;
pub type Topology = topol::Topology;
