use alum::{
    Adaptor, CrossProductAdaptor, DotProductAdaptor, FloatScalarAdaptor, Handle, HasIterators,
    HasTopology, PolyMeshT, VectorLengthAdaptor, VectorNormalizeAdaptor,
};
use three_d::{
    vec3, Camera, CameraAction, CameraControl, Context, CpuMaterial, CpuMesh, Cull, Event, Gm,
    Indices, InnerSpace, InstancedMesh, Instances, Mat4, Mesh, MetricSpace, PhysicalMaterial,
    Positions, Quat, Srgba, Vec3,
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
    sphere.transform(&Mat4::from_scale(vertex_radius)).unwrap();
    let mut cylinder = CpuMesh::cylinder(10);
    cylinder
        .transform(&Mat4::from_nonuniform_scale(1.0, edge_radius, edge_radius))
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

///
/// A control that makes the camera orbit around a target.
///
pub struct CameraMouseControl {
    control: CameraControl,
}

impl CameraMouseControl {
    /// Creates a new orbit control with the given target and minimum and maximum distance to the target.
    pub fn new(target: Vec3, min_distance: f32, max_distance: f32) -> Self {
        Self {
            control: CameraControl {
                right_drag_horizontal: CameraAction::OrbitLeft { target, speed: 0.1 },
                right_drag_vertical: CameraAction::OrbitUp { target, speed: 0.1 },
                left_drag_horizontal: CameraAction::Left { speed: 0.005 },
                left_drag_vertical: CameraAction::Up { speed: 0.005 },
                scroll_vertical: CameraAction::Zoom {
                    min: min_distance,
                    max: max_distance,
                    speed: 0.1,
                    target,
                },
                ..Default::default()
            },
        }
    }

    /// Handles the events. Must be called each frame.
    pub fn handle_events(&mut self, camera: &mut Camera, events: &mut [Event]) -> bool {
        if let CameraAction::Zoom { speed, target, .. } = &mut self.control.scroll_vertical {
            *speed = 0.01 * target.distance(*camera.position()) + 0.001;
        }
        if let CameraAction::OrbitLeft { speed, target } = &mut self.control.right_drag_horizontal {
            *speed = 0.01 * target.distance(*camera.position()) + 0.001;
        }
        if let CameraAction::OrbitUp { speed, target } = &mut self.control.right_drag_vertical {
            *speed = 0.01 * target.distance(*camera.position()) + 0.001;
        }
        if let CameraAction::Left { speed } = &mut self.control.left_drag_horizontal {
            *speed = 0.0005 * camera.target().distance(*camera.position()) + 0.00001;
        }
        if let CameraAction::Up { speed } = &mut self.control.left_drag_vertical {
            *speed = 0.0005 * camera.target().distance(*camera.position()) + 0.00001;
        }
        self.control.handle_events(camera, events)
    }
}
