use alum::{
    Adaptor, CrossProductAdaptor, FloatScalarAdaptor, Handle, HasIterators, HasTopology, PolyMeshT,
    VectorNormalizeAdaptor,
};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, Context, CpuMaterial, CpuMesh, Cull,
    DirectionalLight, FrameOutput, Gm, Indices, InnerSpace, InstancedMesh, Instances, Mat4, Mesh,
    OrbitControl, PhysicalMaterial, Positions, Quat, Srgba, Vec3, Window, WindowSettings,
};

struct MeshAdaptor;

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

type PolygonMesh = PolyMeshT<3, MeshAdaptor>;

fn mesh_view(
    mut mesh: PolygonMesh,
    context: &Context,
    vertex_radius: f32,
    edge_radius: f32,
) -> (
    Gm<Mesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
) {
    let points = mesh.points();
    let points = points.try_borrow().expect("Cannot borrow points");
    let vnormals = mesh.update_vertex_normals_accurate().unwrap();
    let vnormals = vnormals.try_borrow().unwrap();
    let cpumesh = CpuMesh {
        positions: Positions::F32(points.iter().map(|p| vec3(p.x, p.y, p.z)).collect()),
        indices: Indices::U32(
            mesh.triangulated_vertices()
                .flatten()
                .map(|v| v.index())
                .collect(),
        ),
        normals: Some(vnormals.to_vec()),
        ..Default::default()
    };
    let model_material = PhysicalMaterial::new_opaque(
        &context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(200, 200, 200),
            roughness: 0.7,
            metallic: 0.8,
            ..Default::default()
        },
    );
    let mut wireframe_material = PhysicalMaterial::new_opaque(
        &context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(220, 50, 50),
            roughness: 0.7,
            metallic: 0.8,
            ..Default::default()
        },
    );
    wireframe_material.render_states.cull = Cull::Back;
    let mut sphere = CpuMesh::sphere(8);
    sphere.transform(&Mat4::from_scale(vertex_radius)).unwrap();
    let mut cylinder = CpuMesh::cylinder(10);
    cylinder
        .transform(&Mat4::from_nonuniform_scale(1.0, edge_radius, edge_radius))
        .unwrap();
    (
        Gm::new(Mesh::new(&context, &cpumesh), model_material),
        Gm::new(
            InstancedMesh::new(
                &context,
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
                &context,
                &Instances {
                    transformations: mesh
                        .edges()
                        .map(|e| {
                            let h = e.halfedge(false);
                            let mut ev = mesh.calc_halfedge_vector(h, &points);
                            let length = ev.magnitude();
                            ev /= length;
                            let ev = vec3(ev.x, ev.y, ev.z);
                            let start = points[h.tail(&mesh)];
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

fn main() {
    // Window and context.
    let window = Window::new(WindowSettings {
        title: "Viewer".to_string(),
        min_size: (512, 256),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();
    // Setup the camera and the controls and lights.
    let target = vec3(0.5, 0.5, 0.5);
    let scene_radius: f32 = 6.0;
    let mut camera = Camera::new_perspective(
        window.viewport(),
        target + scene_radius * vec3(0.6, 0.3, 1.0).normalize(),
        target,
        vec3(0.0, 1.0, 0.0),
        degrees(45.0),
        0.1,
        1000.0,
    );
    let mut control = OrbitControl::new(*camera.target(), 0.1 * scene_radius, 100.0 * scene_radius);
    let ambient = AmbientLight::new(&context, 0.7, Srgba::WHITE);
    let directional0 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(-1.0, -1.0, -1.0));
    let directional1 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(1.0, 1.0, 1.0));
    // Create the mesh.
    let mesh = PolygonMesh::unit_box().unwrap();
    let (mesh, vertices, edges) = mesh_view(mesh, &context, 0.01, 0.005);
    // Render loop.
    window.render_loop(move |mut frame_input| {
        let mut redraw = frame_input.first_frame;
        redraw |= camera.set_viewport(frame_input.viewport);
        redraw |= control.handle_events(&mut camera, &mut frame_input.events);
        if redraw {
            frame_input
                .screen()
                .clear(ClearState::color_and_depth(0.1, 0.1, 0.1, 1.0, 1.0))
                .render(
                    &camera,
                    mesh.into_iter().chain(&vertices).chain(&edges),
                    &[&ambient, &directional0, &directional1],
                );
        }
        FrameOutput {
            swap_buffers: redraw,
            ..Default::default()
        }
    });
}
