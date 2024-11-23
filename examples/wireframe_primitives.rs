mod common;

use alum::HasTopology;
use common::{mesh_view, PolygonMesh};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, Context, CpuMaterial, CpuMesh, Cull,
    DirectionalLight, FrameOutput, Gm, InnerSpace, InstancedMesh, Instances, Mat4, OrbitControl,
    PhysicalMaterial, Quat, Srgba, Window, WindowSettings,
};

fn wireframe_view(
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
    // Create the meshes.
    let mut tet = PolygonMesh::tetrahedron(1.0).unwrap();
    tet.translate(vec3(-4.0, 0.0, 0.0)).unwrap();
    tet.update_vertex_normals_accurate().unwrap();
    let mut hex = PolygonMesh::hexahedron(1.0).unwrap();
    hex.translate(vec3(-2.0, 0.0, 0.0)).unwrap();
    hex.update_vertex_normals_accurate().unwrap();
    let mut oct = PolygonMesh::octahedron(1.0).unwrap();
    oct.translate(vec3(0.0, 0.0, 0.0)).unwrap();
    oct.update_vertex_normals_accurate().unwrap();
    let mut dod = PolygonMesh::dodecahedron(1.0).unwrap();
    dod.translate(vec3(2.0, 0.0, 0.0)).unwrap();
    dod.update_vertex_normals_accurate().unwrap();
    let mut ico = PolygonMesh::icosahedron(1.0).unwrap();
    ico.translate(vec3(4.0, 0.0, 0.0)).unwrap();
    ico.update_vertex_normals_accurate().unwrap();
    let (tm, (tv, te)) = (
        mesh_view(&tet, &context),
        wireframe_view(&tet, &context, 0.01, 0.005),
    );
    let (hm, (hv, he)) = (
        mesh_view(&hex, &context),
        wireframe_view(&hex, &context, 0.01, 0.005),
    );
    let (om, (ov, oe)) = (
        mesh_view(&oct, &context),
        wireframe_view(&oct, &context, 0.01, 0.005),
    );
    let (dm, (dv, de)) = (
        mesh_view(&dod, &context),
        wireframe_view(&dod, &context, 0.01, 0.005),
    );
    let (im, (iv, ie)) = (
        mesh_view(&ico, &context),
        wireframe_view(&ico, &context, 0.01, 0.005),
    );
    {
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
                        tm.into_iter()
                            .chain(&tv)
                            .chain(&te)
                            .chain(&hm)
                            .chain(&hv)
                            .chain(&he)
                            .chain(&om)
                            .chain(&ov)
                            .chain(&oe)
                            .chain(&dm)
                            .chain(&dv)
                            .chain(&de)
                            .chain(&im)
                            .chain(&iv)
                            .chain(&ie),
                        &[&ambient, &directional0, &directional1],
                    );
            }
            FrameOutput {
                swap_buffers: redraw,
                ..Default::default()
            }
        });
    }
}
