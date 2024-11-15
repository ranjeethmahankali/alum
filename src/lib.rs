/*!
This is a halfedge based polygon mesh library inspired by
[OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/). Hence this
crate has an API that is very similar to OpenMesh.

# Overview

+ A halfedge datastructure is used to represent the topology of a mesh, i.e. the
  connectivity of vertices, edges and faces.

+ The generic polygon mesh type [`PolyMeshT<DIM, A>`] can be used with custom
  geometric types. These include:

  + Vector/Point type used to represent the positions of vertices, normals of
    faces etc.

  + The scalar type used to represent the areas of faces, lengths of edges, and
    angles between edges etc.

  + Number of spacial dimensions the mesh lives in. Often this is 2 or 3.

+ To use the generic [`PolyMeshT<DIM, A>`] with your own geometric types, you
  must provide an implementation of [`Adaptor`] that tells this crate how to
  work with your geometric types. A basic implementation of the [`Adaptor`] is
  needed to use [`PolyMeshT<DIM, A>`]. There are ancilliary adaptor types that
  can be optionally implemented to enable variuos geometric features of
  [`PolyMeshT<DIM, A>`] such as computing areas of faces, normals of faces,
  normals of vertices, dihedral angles etc.

+ Optionally, this crate provides some builtin adaptor implementations, and
  concrete mesh types that can be used without any boilerplate. These builtin
  geometric types use the [`glam`](https://crates.io/crates/glam) crate and can
  be found in the [`alum_glam`] module. `use_glam` feature is required by
  these. These builtin types include:

  + [`PolyMeshF32`](alum_glam::PolyMeshF32) using the
    [`BuiltInAdaptorF32`](alum_glam::BuiltInAdaptorF32). This mesh uses 32 bit
    floating point types to represent the geometry.

  + [`PolyMeshF64`](alum_glam::PolyMeshF64) using the
    [`BuiltInAdaptorF64`](alum_glam::BuiltInAdaptorF64). This mesh uses 64 bit
    floating point types to represent the geometry.
*/

mod check;
mod create;
mod edit;
mod element;
mod error;
mod iterator;
mod macros;
mod math;
mod mesh;
mod obj;
mod property;
mod status;
mod subidv;
mod topol;

#[cfg(feature = "use_glam")]
pub mod alum_glam;

pub use element::{Handle, HasTopology, EH, FH, HH, VH};
pub use error::Error;
pub use mesh::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, PolyMeshT,
    VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use property::{EProperty, FProperty, HProperty, Property, VProperty};
pub use status::Status;
