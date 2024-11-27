# Render | [code](https://github.com/ranjeethmahankali/alum/tree/main/examples/render.rs)

This example shows how to render the triangles and the wireframe of an `alum`
mesh with `three_d`.

```sh
cargo run --example render
```

![Screen capture of a cube](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/render.png)

# Primitives |
[code](https://github.com/ranjeethmahankali/alum/tree/main/examples/primitives.rs)

This example shows how to create primitive polyhedral meshes, and render them
with `three_d`.

```sh
cargo run --example primitives
```

![Screen capture of polyhedra](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/primitives.png)

# Loading an OBJ File | [code](https://github.com/ranjeethmahankali/alum/tree/main/examples/load_obj.rs)

This example shows how to load a mesh from an OBJ file and render it using
`three_d`.

```sh
cargo run --example load_obj
```

![Screen capture of the imported OBJ mesh](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/load_obj.png)

# Property System | [code](https://github.com/ranjeethmahankali/alum/tree/main/examples/property.rs)

This example finds the vertices with extreme coordinate values along x, y, and z
directions. These are considered the 'seed' vertices. The remaining vertices are
then classified according to their closest seed vertex. The distance between the
vertices is measured topological, i.e. the number of edges along the shortest
chain of edges connecting the vertices. The vertices are colored according to
their classification.

This example uses the property system to keep track of the classification and
the colors.

```sh
cargo run --example property
```

![Screen capture of classified mesh vertices](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/property.png)

# Mesh Subdivision | [code](https://github.com/ranjeethmahankali/alum/tree/main/examples/subdiv.rs)

This example demonstrates mesh subdivision using `alum`. It shows three meshes,
subdivided using Catmull-Clark, Loop, and Sqrt-3 subdivision schemes
respectively. You can control the number of subdivisions with up and down
arrows. For performance reasons, it is recommended running this example in
release mode:

```sh
cargo run --release --example subdiv
```

![Screen capture of mesh subdivision](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/subdiv.png)

# Mesh Decimation | [code](https://github.com/ranjeethmahankali/alum/tree/main/examples/decimate.rs)

This example demonstrates mesh decimation using `alum`. You can control the
number of edge collapses using the up arrow (+1), down arrow (-1), page-up
(+10), page-down (-10), and press R to revert to the original mesh. For
performance reasons, it is recommended running this example in release mode:

```sh
cargo run --release --example decimate
```

![Screen capture of mesh decimation](https://github.com/ranjeethmahankali/alum/tree/main/examples/images/decimate.png)
