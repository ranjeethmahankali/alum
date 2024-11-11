mod edit;
mod element;
mod error;
mod iterator;
mod macros;
mod math;
mod mesh;
mod obj;
mod primitive;
mod property;
mod status;
mod topol;

#[cfg(feature = "use_glam")]
pub mod alum_glam;

pub use element::{Handle, EH, FH, HH, VH};
pub use error::Error;
pub use mesh::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, PolyMeshT,
    VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use property::{EProperty, FProperty, HProperty, Property, VProperty};
pub use status::Status;
pub use topol::Topology;
