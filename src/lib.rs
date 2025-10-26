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

+ This crate provides some default adaptor implementations, and concrete mesh
  types that can be used without any boilerplate. These use trivial vector and
  implementation.

  + [`PolyMeshF32`] uses 32 bit floating point types to represent the geometry.

  + [`PolyMeshF64`] uses 64 bit floating point types to represent the geometry.

# Halfedge Mesh Representation

```text
   ^|          ^|          ^|          ^|
   ||          ||          ||          ||
   ||          ||          ||          ||
   ||          ||          ||          ||
   ||          ||          ||          ||
   |v          |v          |v          |v
   ** -------> ** -------> ** -------> **
   ** <------- ** <------- ** <------- **
   ^|          ^|          ^|          ^|
   ||          ||          ||          ||
   ||          ||          ||          ||
   ||          ||          ||          ||
   ||          ||          ||          ||
   |v          |v          |v          |v
```

Each is represented as a [`pair of halfedges`](EH::halfedges), that point in
opposite directions. These two halfedges are considered the
[`opposite`](HH::opposite) of each other. The vertex a halfedge points toward is
its [`head`](HH::head) and the vertex it points away from is its
[`tail`](HH::tail). So for a pair of opposite halfedges belonging to an edge,
the head of one of them is the tail of the other and vice-versa. A circular
chain of halfedges where each consecutive halfedge points away from the head of
the previous halfedge form a `loop`. A halfedge must always exist in a loop. The
halfedges in a loop are linked to each other as [`next`](HH::next) and
[`previous`](HH::prev) halfedges. A loop can optionally contain a face,
otherwise it's considered a boundary loop. When a loop is not a boundary loop,
every halfedge in that loop is linked to the [`face`](HH::face), otherwise the
halfedge is considered a [boundary](HH::is_boundary). If either halfedge in an
edge is boundary, then the edge is also considered [boundary](EH::is_boundary).

Every vertex in the mesh is optionally linked to a
[`halfedge`](VH::halfedge). If a vertex is not linked to a halfedge, it is
considered isolated. When a vertex is linked to a halfedge, all other topology
incident on that vertex can be reached by traversing from that halfedge. If a
vertex is on the boundary of a mesh, it must be linked to a boundary
halfedge. This assumption is used to check if a vertex [is on the
boundary](VH::is_boundary).

Every face is linked to exactly one [`halfedge`](FH::halfedge) from it's
loop. All other halfedges in that loop, can be reached by traversing the
[`next`](HH::next) and [`previous`](HH::prev) links of the halfedge.

For more on how to work with halfedge mesh datastructures, even though the
vocabulary used by OpenMesh is slightly different it is identical to the one
used by this crate and can be super helpful:

- <https://www.graphics.rwth-aachen.de/media/openmesh_static/Documentations/OpenMesh-6.0-Documentation/a00016.html>
- <https://www.graphics.rwth-aachen.de/media/openmesh_static/Documentations/OpenMesh-6.2-Documentation/a00032.html>

In this crate, these elements are addressed using handles: [`VH`] for vertex,
[`HH`] for halfedge, [`EH`] for edge, and [`FH`] for face. The handles are just
convenient wrappers around the indices of these elements. So when traversing
topology, the container with the topological information (e.g. the mesh) must be
supplied as shown in below example. While this is true for most topological
functions, note that some of them don't require the mesh to be passed in. For
example, the [opposite](HH::opposite) of a halfedge, and the
[halfedge](EH::halfedge) of an edge, are implicitly computed based on it's
indices. This is a consequence of how the halfedges and edges are stored.

```rust
use alum::{PolyMeshF32, HasTopology, HasIterators, Handle};

let mesh = PolyMeshF32::unit_box().expect("Cannot create a box");

// Find halfedge pointing from vertex 0 to vertex 1.
let h = mesh.find_halfedge(0.into(), 1.into()).expect("Cannot find halfedge");

// To traverse to the tail vertex or the head vertex of this halfedge, we have to
// provide the mesh because `h` is just a thin wrapper around the index.
let v0 = h.tail(&mesh);
let v1 = h.head(&mesh);
assert_eq!(v0.index(), 0);
assert_eq!(v1.index(), 1);

// The traversal functions can be chained to 'walk` the mesh.
let f = h.next(&mesh).next(&mesh).opposite().face(&mesh);
assert!(f.is_some()); // No boundary halfedges in a box. So this must link to a face.
```
 */

mod check;
mod create;
mod edit;
mod element;
mod error;
mod history;
mod iterator;
mod macros;
mod math;
mod mesh;
mod obj;
mod property;
mod queue;
mod status;
mod topol;

#[cfg(feature = "subdiv")]
pub mod subdiv;

#[cfg(feature = "decimate")]
pub mod decimate;

#[cfg(feature = "decimate")]
pub use decimate::{
    Decimater, HasDecimation, edge_length::EdgeLengthDecimater, quadric::Quadric,
    quadric::QuadricDecimater, quadric::QuadricType,
};

pub use edit::EditableTopology;
pub use element::{EH, FH, HH, Handle, HandleRange, VH};
pub use error::Error;
pub use history::{
    EPropHistory, FPropHistory, HPropHistory, PropHistory, TopolHistory, VPropHistory,
};
pub use iterator::HasIterators;
pub use mesh::{
    Adaptor, CrossProductAdaptor, DVec3, DotProductAdaptor, FloatScalarAdaptor, PolyMeshF32,
    PolyMeshF64, PolyMeshT, Vec3, VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use property::{
    EPropBuf, EProperty, FPropBuf, FProperty, HPropBuf, HProperty, PropBuf, Property, VPropBuf,
    VProperty,
};
pub use queue::Queue;
pub use status::Status;
pub use topol::HasTopology;
