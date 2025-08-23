# Alum: A Halfedge based Polygon Mesh Library

[![CI](https://github.com/ranjeethmahankali/alum/actions/workflows/ci.yml/badge.svg)](https://github.com/ranjeethmahankali/alum/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.rs/alum/latest/alum/)
[![Crate](https://img.shields.io/crates/v/alum)](https://crates.io/crates/alum)
[![Examples with three_d](https://img.shields.io/badge/three__d-examples-purple)](https://github.com/ranjeethmahankali/alum/tree/main/examples)

![Standford Bunny](assets/bunny.png)

This library is inspired by
[OpenMesh](https://www.graphics.rwth-aachen.de/software/openmesh/), hence has an
API very similar to that of OpenMesh. I love using OpenMesh in C++, and wrote
this library because I couldn't find an equivalent in Rust. Huge thanks to
OpenMesh and it's maintainers for the inspiration!

For now, the features of this library exist to serve my other projects. It has
full parity with the `OpenMesh/Core` module, and has partial parity with the
`OpenMesh/Tools` module, and new features will be added as required by my other
projects.

## Installation

This can be added to a Rust project as a dependency from
[crates.io](https://crates.io/crates/alum) with:

```sh
cargo add alum
```

Or by adding the following to your `Cargo.toml`:

```toml
[dependencies]
alum = "0.6.1"
```

## Usage and Features

The library includes default implementations of mesh adaptors, and corresponding
mesh types with single and double floating point precision. If you're not happy
with these default adaptors, you can write your own adaptor that tells the
library how to work with the vector and scalar types of your choice. Read the
[documentation](https://docs.rs/alum/latest/alum/) to learn more about
this. These
[examples](https://github.com/ranjeethmahankali/alum/tree/main/examples)
demonstrate writing custom adaptors and rendering and using various features of
this crate together with [`three_d`](https://github.com/asny/three-d)

### Subdivision

Mesh subdivision is available under the "subdiv" feature. Catmull-Clark, Loop,
and Sqrt3 subdivision schemes are implemented. [This
example](https://github.com/ranjeethmahankali/alum/tree/main/examples)
demonstrates all three decimation schemes.

### Decimation

Decimation is implemented with an architecture very similar to that of OpenMesh,
and is available under the "decimate" feature. You can write your own
implementation of the `Decimater` trait, to decimate the mesh with a custom
heuristic. The library does come with some commonly useful decimater
implementations. More maybe added later. [This
example](https://github.com/ranjeethmahankali/alum/tree/main/examples)
demonstrates the decimation of a mesh using the probabilistic quadric minimizing decimater.

### Property System

This library also comes with a property system just like the one in `OpenMesh`,
with some small improvements and differences. The properties are always
synchronized with the mesh elements, through additions, deletions and garbage
collections which result in reordering of elements. Unlike properties in
`OpenMesh`, the properties here use interior mutability with `RefCell<T>` and
enforce runtime borrow checking rules. Read the
[documentation](https://docs.rs/alum/latest/alum/) to learn more.
