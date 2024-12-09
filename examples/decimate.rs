mod common;

use alum::{
    Decimater, Handle, HasDecimation, HasIterators, HasTopology, QuadricDecimater, QuadricType,
};
use common::{CameraMouseControl, MeshAdaptor, PolygonMesh};
use std::{path::PathBuf, time::Instant};
use three_d::{
    degrees, vec3, AmbientLight, Camera, ClearState, Context, CpuMaterial, CpuMesh, Cull,
    DirectionalLight, Event, FrameOutput, Gm, Indices, InnerSpace, InstancedMesh, Instances, Key,
    Mat4, Mesh, PhysicalMaterial, Positions, Quat, Srgba, Viewport, Window, WindowSettings,
};

pub struct RecordingDecimater {
    inner: QuadricDecimater<MeshAdaptor>,
    history: Vec<PolygonMesh>,
}

impl RecordingDecimater {
    pub fn new(mesh: &PolygonMesh) -> Self {
        RecordingDecimater {
            inner: QuadricDecimater::new(mesh, QuadricType::ProbabilisticTriangle)
                .expect("Cannot create quadric decimater"),
            history: Vec::new(),
        }
    }

    pub fn history(&self) -> &[PolygonMesh] {
        &self.history
    }
}

impl Decimater<PolygonMesh> for RecordingDecimater {
    type Cost = f32;

    fn collapse_cost(&self, mesh: &PolygonMesh, h: alum::HH) -> Option<Self::Cost> {
        self.inner.collapse_cost(mesh, h)
    }

    fn before_collapse(&mut self, mesh: &PolygonMesh, h: alum::HH) -> Result<(), alum::Error> {
        self.inner.before_collapse(mesh, h)?;
        Ok(())
    }

    fn after_collapse(&mut self, mesh: &PolygonMesh, v: alum::VH) -> Result<(), alum::Error> {
        self.inner.after_collapse(mesh, v)?;
        let mut copy = mesh.try_clone()?;
        copy.garbage_collection()?;
        copy.check_topology()?;
        self.history.push(copy);
        Ok(())
    }
}

fn bunny_mesh() -> PolygonMesh {
    let mesh = PolygonMesh::load_obj(PathBuf::from("/home/rnjth94/dev/alum/assets/bunny.obj"))
        .expect("Cannot load obj");
    {
        let mut points = mesh.points();
        let mut points = points.try_borrow_mut().expect("Cannot borrow points");
        for p in points.iter_mut() {
            *p *= 10.; // Scale the mesh.
        }
    }
    mesh
}

fn visualize_mesh(mesh: &PolygonMesh, context: &Context) -> MeshViews {
    const CYL_RADIUS: f32 = 0.0005;
    const SPH_RADIUS: f32 = 0.001;
    // Create a CPU-side mesh consisting of a single colored triangle
    let (model, etransforms, vtransforms) = {
        let points = mesh.points();
        let points = points.try_borrow().expect("Cannot borrow points");
        let vnormals = mesh.vertex_normals().expect("Cannot borrow vertex normals");
        let vnormals = vnormals.try_borrow().expect("Cannot borrow vertex normals");
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
            normals: Some(vnormals.iter().map(|n| vec3(n.x, n.y, n.z)).collect()),
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
        (
            Gm::new(Mesh::new(context, &cpumesh), model_material),
            Instances {
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
            Instances {
                transformations: points
                    .iter()
                    .map(|pos| Mat4::from_translation(vec3(pos.x, pos.y, pos.z)))
                    .collect(),
                ..Default::default()
            },
        )
    };
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
    sphere.transform(&Mat4::from_scale(SPH_RADIUS)).unwrap();
    let vertices = Gm::new(
        InstancedMesh::new(context, &vtransforms, &sphere),
        wireframe_material.clone(),
    );
    let mut cylinder = CpuMesh::cylinder(10);
    cylinder
        .transform(&Mat4::from_nonuniform_scale(1.0, CYL_RADIUS, CYL_RADIUS))
        .unwrap();
    let edges = Gm::new(
        InstancedMesh::new(context, &etransforms, &cylinder),
        wireframe_material,
    );
    MeshViews(model, vertices, edges)
}

struct MeshViews(
    Gm<Mesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
    Gm<InstancedMesh, PhysicalMaterial>,
);

pub fn main() {
    let window = Window::new(WindowSettings {
        title: "Viewer".to_string(),
        min_size: (512, 256),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();
    let mut mesh = {
        let mut mesh = bunny_mesh();
        mesh.update_face_normals()
            .expect("Cannot update face normals");
        mesh.update_vertex_normals_fast()
            .expect("Cannot update vertex normals");
        mesh.check_topology().expect("Topological errors found");
        mesh
    };
    let mut decimater = RecordingDecimater::new(&mesh);
    const NUM_COLLAPSES: usize = 3000;
    let before = Instant::now();
    mesh.decimate(&mut decimater, NUM_COLLAPSES)
        .expect("Cannot decimate");
    println!(
        "Decimation took {}ms",
        (Instant::now() - before).as_millis()
    );
    let history = decimater.history();
    let target = vec3(0.0, 1.0, 0.0);
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
    let views: Vec<MeshViews> = history
        .iter()
        .map(|mesh| visualize_mesh(mesh, &context))
        .collect();
    // let (_, rvs, res) = visualize_mesh(&refbox, &context);
    let ambient = AmbientLight::new(&context, 0.7, Srgba::WHITE);
    let directional0 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(-1.0, -1.0, -1.0));
    let directional1 = DirectionalLight::new(&context, 2.0, Srgba::WHITE, &vec3(1.0, 1.0, 1.0));
    let mut index = 0usize;
    let mut gui = three_d::GUI::new(&context);
    // render loop
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
                        ui.label("Use the Up/Down arrows to control the number of edge collapses.");
                        ui.label("Use PageUp / PageDown to control edge collapses in steps of 10.");
                        ui.label("Press R to revert back to the original mesh.");
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
            let next = usize::clamp(
                match event {
                    Event::KeyPress { kind, .. } => match kind {
                        Key::ArrowUp => index + 1,
                        Key::ArrowDown => index.saturating_sub(1),
                        Key::R => 0,
                        Key::PageUp => index + 10,
                        Key::PageDown => index.saturating_sub(10),
                        _ => index,
                    },
                    _ => index,
                },
                0,
                views.len() - 1,
            );
            if next != index {
                redraw = true;
                index = next;
            }
        }
        let MeshViews(model, vertices, edges) = &views[index];
        if redraw {
            frame_input
                .screen()
                .clear(ClearState::color_and_depth(0.1, 0.1, 0.1, 1.0, 1.0))
                .render(
                    &camera,
                    model.into_iter().chain(vertices).chain(edges),
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
