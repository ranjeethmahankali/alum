use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, PolyMeshT, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
use three_d::{
    Context, CpuMaterial, CpuMesh, Cull, Gm, Indices, InnerSpace, InstancedMesh, Instances, Mat4,
    Mesh, PhysicalMaterial, Positions, Quat, Srgba, Vec3, vec3,
};

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

impl VectorLengthAdaptor<3> for MeshAdaptor {
    fn vector_length(v: Self::Vector) -> Self::Scalar {
        v.magnitude()
    }

    fn vector_length_squared(v: Self::Vector) -> Self::Scalar {
        v.magnitude2()
    }
}

impl DotProductAdaptor<3> for MeshAdaptor {
    fn dot_product(a: Self::Vector, b: Self::Vector) -> Self::Scalar {
        a.dot(b)
    }
}

pub type PolygonMesh = PolyMeshT<3, MeshAdaptor>;

#[allow(dead_code)]
pub fn mesh_view(mesh: &PolygonMesh, context: &Context) -> Gm<Mesh, PhysicalMaterial> {
    let points = mesh.points();
    let points = points.try_borrow().expect("Cannot borrow points");
    let vnormals = mesh.vertex_normals().unwrap();
    let vnormals = vnormals.try_borrow().unwrap();
    let fstatus = mesh.face_status_prop();
    let fstatus = fstatus.try_borrow().unwrap();
    let cpumesh = CpuMesh {
        positions: Positions::F32(points.iter().map(|p| vec3(p.x, p.y, p.z)).collect()),
        indices: Indices::U32(
            mesh.triangulated_vertices(&fstatus)
                .flatten()
                .map(|v| v.index())
                .collect(),
        ),
        normals: Some(vnormals.to_vec()),
        ..Default::default()
    };
    let model_material = PhysicalMaterial::new_opaque(
        context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(200, 200, 200),
            roughness: 0.7,
            metallic: 0.8,
            ..Default::default()
        },
    );
    Gm::new(Mesh::new(context, &cpumesh), model_material)
}

#[allow(dead_code)]
pub fn wireframe_view(
    mesh: &PolygonMesh,
    context: &Context,
    vertex_radius: f32,
    edge_radius: f32,
) -> (
    Gm<InstancedMesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
) {
    let points = mesh.points();
    let points = points.try_borrow().expect("Cannot borrow points");
    let mut wireframe_material = PhysicalMaterial::new_opaque(
        context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(220, 50, 50),
            roughness: 0.7,
            metallic: 0.8,
            ..Default::default()
        },
    );
    wireframe_material.render_states.cull = Cull::Back;
    let mut sphere = CpuMesh::sphere(8);
    sphere.transform(Mat4::from_scale(vertex_radius)).unwrap();
    let mut cylinder = CpuMesh::cylinder(10);
    cylinder
        .transform(Mat4::from_nonuniform_scale(1.0, edge_radius, edge_radius))
        .unwrap();
    (
        Gm::new(
            InstancedMesh::new(
                context,
                &Instances {
                    transformations: points
                        .iter()
                        .map(|pos| Mat4::from_translation(vec3(pos.x, pos.y, pos.z)))
                        .collect(),
                    ..Default::default()
                },
                &sphere,
            ),
            wireframe_material.clone(),
        ),
        Gm::new(
            InstancedMesh::new(
                context,
                &Instances {
                    transformations: mesh
                        .edges()
                        .map(|e| {
                            let h = e.halfedge(false);
                            let mut ev = mesh.calc_halfedge_vector(h, &points);
                            let length = ev.magnitude();
                            ev /= length;
                            let ev = vec3(ev.x, ev.y, ev.z);
                            let start = points[h.tail(mesh)];
                            let start = vec3(start.x, start.y, start.z);
                            Mat4::from_translation(start)
                                * Into::<Mat4>::into(Quat::from_arc(vec3(1.0, 0., 0.0), ev, None))
                                * Mat4::from_nonuniform_scale(length, 1., 1.)
                        })
                        .collect(),
                    ..Default::default()
                },
                &cylinder,
            ),
            wireframe_material,
        ),
    )
}
