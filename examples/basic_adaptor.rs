use alum::{Adaptor, CrossProductAdaptor, FloatScalarAdaptor, PolyMeshT, VectorNormalizeAdaptor};
use three_d::{InnerSpace, Vec3};

pub struct MeshAdaptor;

impl Adaptor<3> for MeshAdaptor {
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

impl FloatScalarAdaptor<3> for MeshAdaptor {
    fn scalarf32(val: f32) -> Self::Scalar {
        val
    }

    fn scalarf64(val: f64) -> Self::Scalar {
        val as f32
    }
}

impl CrossProductAdaptor for MeshAdaptor {
    fn cross_product(a: Self::Vector, b: Self::Vector) -> Self::Vector {
        a.cross(b)
    }
}

impl VectorNormalizeAdaptor<3> for MeshAdaptor {
    fn normalized_vec(v: Self::Vector) -> Self::Vector {
        v.normalize()
    }
}

pub type PolygonMesh = PolyMeshT<3, MeshAdaptor>;

// fn make_box() -> PolygonMesh {
//     // Make a box.
//     let mut mesh = PolygonMesh::with_capacity(8, 12, 6);
//     // Add several vertices at once.
//     let mut vs: Vec<_> = mesh
//         .add_vertices(&[
//             vec3(0.0, 0.0, 0.0),
//             vec3(1.0, 0.0, 0.0),
//             vec3(1.0, 1.0, 0.0),
//             vec3(0.0, 1.0, 0.0),
//             vec3(0.0, 0.0, 1.0),
//             vec3(1.0, 0.0, 1.0),
//         ])
//         .unwrap()
//         .collect();
//     // Add vertices one at a time.
//     let v6 = mesh.add_vertex(vec3(1.0, 1.0, 1.0)).unwrap();
//     let v7 = mesh.add_vertex(vec3(0.0, 1.0, 1.0)).unwrap();
//     vs.push(v6);
//     vs.push(v7);
//     // Add quad faces from 4 vertices.
//     mesh.add_quad_face(vs[0], vs[3], vs[2], vs[1]).unwrap();
//     mesh.add_quad_face(vs[0], vs[1], vs[5], vs[4]).unwrap();
//     mesh.add_quad_face(vs[1], vs[2], vs[6], vs[5]).unwrap();
//     // Add polygon faces from a slice of vertices.
//     mesh.add_face(&[vs[2], vs[3], vs[7], vs[6]]).unwrap();
//     mesh.add_face(&[vs[3], vs[0], vs[4], vs[7]]).unwrap();
//     mesh.add_face(&[vs[4], vs[5], vs[6], vs[7]]).unwrap();
//     mesh
// }
