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

/*!
# Working with mutable iterators

Iterators such as [`Self::loop_ccw_iter`] and other `*_iter` iterators borrow
the mesh immutably. These are useful when iterating over the elements of the
mesh without needing to modify the mesh. While the mesh is borrowed by the
iterators, the borrow checker will not let you borrow the mesh again mutably,
for as long as the iterator is alive. This is a problem if you're trying to
modify the mesh while iterating over it's elements. In such scenarios, mutable
iterators are useful. These functions end with `*_iter_mut`.

The `*_iter_mut` functions borrow the mesh mutably. Similar to immutable
iterators, while the mutable iterator is alive, the borrow checker won't let you
mutably borrow the mesh again. But, unlike the immutable iterators, the mutable
iterators don't yield just the elements of the mesh. Instead the mutable
iterators yield a tuple containing a mutable reference to the mesh, and the mesh
element. You can modify the mesh using the mutable reference yielded by the
mutable iterator. So essentially, the borrow checker is happy because the
iterator borrows the mesh mutably and to modify the mesh, you borrow the mesh
from the iterator. The borrows are chained instead of being simultaenous. This
keeps the borrow checker happy, and ensures safety to some extent. Below is some
example code that iterates over the vertices of face with index 2, and modifies
their positions.

```rust
use alum::{alum_glam::PolyMeshF32, FH, HasIterators};

let mut boxmesh = PolyMeshF32::unit_box().expect("Cannot create a box");
assert_eq!(1.0, boxmesh.try_calc_volume().expect("Cannot compute volume"));
// Modify the mesh while using a mutable iterator - pull points closer to origin.
let f: FH = 2.into();
for (mesh, v) in boxmesh.fv_ccw_iter_mut(f) {
    // Inside this loop, while the iterator is alive, we cannot borrow `boxmesh`
    // because the iterator already borrowed `boxmesh` mutably. Instead we will use
    // the `mesh` mutable reference yielded by the iterator along with the halfedge.
    let mut pos = mesh.point(v).expect("Cannot read position");
    pos *= 0.75;
    mesh.set_point(v, pos);
}
// The iterator is not alive anymore, so we can borrow `boxmesh` again.
boxmesh.check_topology().expect("Topological errors found");
// Volume is smaller than one because we pulled the vertices closer to the origin.
assert!(1.0 > boxmesh.try_calc_volume().expect("Cannot compute volume"));
```

Despite enforcing the borrow checker rules, the mutable iterators can lead to
problems when used incorrectly. Modifying the topology of the mesh while using
mutable iterators is NOT advised, as this can lead to topological errors. Only
do this if you're know what you're doing. This is akin to mutably iterating over
a linked list while modifying the links between the elements of the linked
list. You can easily create cycles, infinite loops and other problems if you're
not careful.
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
mod subdiv;
mod topol;

#[cfg(feature = "use_glam")]
pub mod alum_glam;

pub use edit::EditableTopology;
pub use element::{Handle, EH, FH, HH, VH};
pub use error::Error;
pub use iterator::HasIterators;
pub use mesh::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, PolyMeshT,
    VectorAngleAdaptor, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
pub use property::{EProperty, FProperty, HProperty, Property, VProperty};
pub use status::Status;
pub use topol::HasTopology;
