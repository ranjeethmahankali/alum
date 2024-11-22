use alum::{Adaptor, HasIterators, HasTopology, PolyMeshT};
use three_d::{vec3, Vec3};

pub struct MyAdaptor;

impl Adaptor<3> for MyAdaptor {
    type Vector = Vec3;

    type Scalar = f32;

    fn vector(coords: [Self::Scalar; 3]) -> Self::Vector {
        three_d::vec3(coords[0], coords[1], coords[2])
    }

    fn zero_vector() -> Self::Vector {
        three_d::vec3(0.0, 0.0, 0.0)
    }

    fn vector_coord(v: &Self::Vector, i: usize) -> Self::Scalar {
        v[i]
    }
}

pub type MyMesh = PolyMeshT<3, MyAdaptor>;

fn main() {
    // Make a Box.
    let mut mesh = MyMesh::with_capacity(8, 12, 6);
    let vs: Vec<_> = mesh
        .add_vertices(&[
            vec3(0.0, 0.0, 0.0),
            vec3(1.0, 0.0, 0.0),
            vec3(1.0, 1.0, 0.0),
            vec3(0.0, 1.0, 0.0),
            vec3(0.0, 0.0, 1.0),
            vec3(1.0, 0.0, 1.0),
            vec3(1.0, 1.0, 1.0),
            vec3(0.0, 1.0, 1.0),
        ])
        .unwrap()
        .collect();
    mesh.add_quad_face(vs[0], vs[3], vs[2], vs[1]).unwrap();
    mesh.add_quad_face(vs[0], vs[1], vs[5], vs[4]).unwrap();
    mesh.add_quad_face(vs[1], vs[2], vs[6], vs[5]).unwrap();
    mesh.add_quad_face(vs[2], vs[3], vs[7], vs[6]).unwrap();
    mesh.add_quad_face(vs[3], vs[0], vs[4], vs[7]).unwrap();
    mesh.add_quad_face(vs[4], vs[5], vs[6], vs[7]).unwrap();
    // Print stats.
    println!(
        "Mesh has {} vertices, {} edges, and {} faces.",
        mesh.num_vertices(),
        mesh.num_edges(),
        mesh.num_faces()
    );
    // Print vertex connectivity.
    let points = mesh.points();
    let points = points.try_borrow().unwrap();
    for (v, pos) in mesh.vertices().zip(points.iter()) {
        println!(
            "Vertex at ({}, {}, {}) is connected to vertices at:",
            pos.x, pos.y, pos.z
        );
        for nv in mesh.vv_ccw_iter(v) {
            let pos = points[nv];
            println!("\t- ({}, {}, {})", pos.x, pos.y, pos.z);
        }
    }
}
