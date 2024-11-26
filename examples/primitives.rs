mod common;

use common::{mesh_view, wireframe_view, CameraMouseControl, PolygonMesh};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, DirectionalLight, FrameOutput, InnerSpace,
    Srgba, Window, WindowSettings,
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
    let mut control =
        CameraMouseControl::new(*camera.target(), 0.1 * scene_radius, 100.0 * scene_radius);
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
