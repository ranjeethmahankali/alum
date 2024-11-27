use std::path::PathBuf;

use common::{mesh_view, wireframe_view, CameraMouseControl, PolygonMesh};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, Context, DirectionalLight, Event, FrameOutput,
    Gm, InnerSpace, InstancedMesh, Key, Mesh, Object, PhysicalMaterial, Srgba, Vec3, Viewport,
    Window, WindowSettings,
};

mod common;

fn bounds(mesh: &PolygonMesh) -> (Vec3, Vec3) {
    let points = mesh.points();
    let points = points.try_borrow().unwrap();
    points.iter().fold(
        (
            vec3(f32::MAX, f32::MAX, f32::MAX),
            vec3(f32::MIN, f32::MIN, f32::MIN),
        ),
        |(min, max), p| {
            (
                vec3(
                    f32::min(min.x, p.x),
                    f32::min(min.y, p.y),
                    f32::min(min.z, p.z),
                ),
                vec3(
                    f32::max(max.x, p.x),
                    f32::max(max.y, p.y),
                    f32::max(max.z, p.z),
                ),
            )
        },
    )
}

fn bunny() -> PolygonMesh {
    let path = {
        let mut path = PathBuf::new();
        path.push(env!("CARGO_MANIFEST_DIR"));
        path.push("assets");
        path.push("bunny.obj");
        path
    };
    let mut mesh = PolygonMesh::load_obj(path).unwrap();
    {
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut().unwrap();
        for pos in points.iter_mut() {
            *pos *= 10.0;
        }
    }
    mesh.update_vertex_normals_accurate().unwrap();
    mesh
}

fn translate(mesh: PolygonMesh, v: Vec3) -> PolygonMesh {
    let mut points = mesh.points();
    let mut points = points.try_borrow_mut().unwrap();
    for pos in points.iter_mut() {
        *pos += v;
    }
    mesh
}

fn get_views(
    mesh: &PolygonMesh,
    context: &Context,
) -> (
    Gm<Mesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
) {
    let mview = mesh_view(mesh, context);
    let (verts, edges) = wireframe_view(mesh, context, 0.001, 0.0005);
    (mview, verts, edges)
}

enum SubdivType {
    CatmullClark,
    Loop,
    Sqrt3(bool),
}

struct MeshData {
    mesh: PolygonMesh,
    stype: SubdivType,
    triangles: Gm<Mesh, PhysicalMaterial>,
    vertices: Gm<InstancedMesh, PhysicalMaterial>,
    edges: Gm<InstancedMesh, PhysicalMaterial>,
}

impl MeshData {
    fn new(mesh: PolygonMesh, stype: SubdivType, context: &Context) -> Self {
        let triangles = mesh_view(&mesh, context);
        let (vertices, edges) = wireframe_view(&mesh, context, 0.001, 0.0005);
        Self {
            mesh,
            stype,
            triangles,
            vertices,
            edges,
        }
    }

    fn into_iter(&self) -> impl Iterator<Item = &dyn Object> {
        self.triangles
            .into_iter()
            .chain(&self.vertices)
            .chain(&self.edges)
    }

    fn chain<'a>(
        &'a self,
        iter: impl Iterator<Item = &'a dyn Object>,
    ) -> impl Iterator<Item = &'a dyn Object> {
        iter.chain(&self.triangles)
            .chain(&self.vertices)
            .chain(&self.edges)
    }

    fn subdivide(&self, context: &Context) -> Self {
        let mut mesh = self.mesh.clone();
        let stype = match self.stype {
            SubdivType::CatmullClark => {
                mesh.subdivide_catmull_clark(1, true).unwrap();
                SubdivType::CatmullClark
            }
            SubdivType::Loop => {
                mesh.subdivide_loop(1, true).unwrap();
                SubdivType::Loop
            }
            SubdivType::Sqrt3(phase) => SubdivType::Sqrt3(!mesh.subdivide_sqrt3(1, phase).unwrap()),
        };
        Self::new(mesh, stype, context)
    }
}

struct MeshTriplet(MeshData, MeshData, MeshData);

impl MeshTriplet {
    fn subdivide(&self, context: &Context) -> Self {
        MeshTriplet(
            self.0.subdivide(context),
            self.1.subdivide(context),
            self.2.subdivide(context),
        )
    }

    fn into_iter(&self) -> impl Iterator<Item = &dyn Object> {
        self.0
            .into_iter()
            .chain(self.1.into_iter())
            .chain(self.2.into_iter())
    }
}

fn main() {
    let window = Window::new(WindowSettings {
        title: "Example: Catmull-Clark, Loop, and Sqrt-3 Subdivision of a Triangle Mesh"
            .to_string(),
        min_size: (512, 256),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();
    let (mut views, target) = {
        let mesh = bunny();
        let (min, max) = bounds(&mesh);
        let target = 0.5 * (min + max);
        const SHIFT: f32 = 2.0;
        (
            vec![MeshTriplet(
                MeshData::new(
                    translate(mesh.clone(), vec3(-SHIFT, 0.0, 0.0)),
                    SubdivType::CatmullClark,
                    &context,
                ),
                MeshData::new(mesh.clone(), SubdivType::Loop, &context),
                MeshData::new(
                    translate(mesh, vec3(SHIFT, 0.0, 0.0)),
                    SubdivType::Sqrt3(false),
                    &context,
                ),
            )],
            target,
        )
    };
    // Setup the camera and the controls and lights.
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
    let mut gui = three_d::GUI::new(&context);
    let mut num_iter = 0usize;
    window.render_loop(move |mut frame_input| {
        let mut panel_height = 0.0;
        gui.update(
            &mut frame_input.events,
            frame_input.accumulated_time,
            frame_input.viewport,
            frame_input.device_pixel_ratio,
            |gui_context| {
                use three_d::egui::*;
                TopBottomPanel::top("top_panel").show(gui_context, |ui| {
                    ui.vertical_centered(|ui| {
                        ui.label("Use the Up/Down arrows to control the number of iterations");
                    });
                });
                panel_height = gui_context.used_rect().height();
            },
        );
        let panel_height = panel_height * frame_input.device_pixel_ratio;
        let viewport = Viewport {
            x: 0,
            y: panel_height as i32,
            width: frame_input.viewport.width,
            height: frame_input.viewport.height - panel_height as u32,
        };
        camera.set_viewport(viewport);
        // First render the view.
        let mut redraw = frame_input.first_frame;
        redraw |= camera.set_viewport(frame_input.viewport);
        redraw |= control.handle_events(&mut camera, &mut frame_input.events);
        for event in &frame_input.events {
            let next = match event {
                Event::KeyPress { kind, .. } => match kind {
                    Key::ArrowDown => num_iter + 1,
                    Key::ArrowUp => num_iter.saturating_sub(1),
                    _ => num_iter,
                },
                _ => num_iter,
            };
            if next != num_iter {
                redraw = true;
                num_iter = next;
            }
            while num_iter >= views.len() {
                views.push(views.last().unwrap().subdivide(&context));
            }
        }
        let view = &views[num_iter];
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
        // Then render the gui on top of the view.
        frame_input.screen().write(|| gui.render()).unwrap();
        FrameOutput {
            swap_buffers: redraw,
            ..Default::default()
        }
    });
}
