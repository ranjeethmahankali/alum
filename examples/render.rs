mod common;

use common::{mesh_view, PolygonMesh};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, DirectionalLight, FrameOutput, InnerSpace,
    OrbitControl, Srgba, Window, WindowSettings,
};

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
    let mut mesh = PolygonMesh::unit_box().unwrap();
    mesh.update_vertex_normals_accurate().unwrap();
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
