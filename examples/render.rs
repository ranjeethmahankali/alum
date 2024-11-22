mod adaptor;

use adaptor::{make_box, PolygonMesh};
use alum::{Handle, HasIterators};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, ColorMaterial, Context, CpuMaterial, CpuMesh,
    DirectionalLight, FrameOutput, Gm, Indices, InnerSpace, Mesh, OrbitControl, Positions, Srgba,
    Window, WindowSettings,
};

fn mesh_view(mesh: &PolygonMesh, context: &Context) -> Gm<Mesh, ColorMaterial> {
    let points = mesh.points();
    let points = points.try_borrow().expect("Cannot borrow points");
    let cpumesh = CpuMesh {
        positions: Positions::F32(points.iter().map(|p| vec3(p.x, p.y, p.z)).collect()),
        indices: Indices::U32(
            mesh.triangulated_vertices()
                .flatten()
                .map(|v| v.index())
                .collect(),
        ),
        ..Default::default()
    };
    let model_material = ColorMaterial::new_opaque(
        &context,
        &CpuMaterial {
            albedo: Srgba::new_opaque(200, 200, 200),
            ..Default::default()
        },
    );
    Gm::new(Mesh::new(&context, &cpumesh), model_material)
}

#[allow(dead_code)]
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
    let mesh = make_box();
    let view = mesh_view(&mesh, &context);
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
                    view.into_iter(),
                    &[&ambient, &directional0, &directional1],
                );
        }
        FrameOutput {
            swap_buffers: redraw,
            ..Default::default()
        }
    });
}
