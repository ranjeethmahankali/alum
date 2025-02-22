mod common;

use alum::{Handle, HasIterators, HasTopology, VH, VProperty};
use common::{CameraMouseControl, PolygonMesh};
use core::f32;
use std::{collections::VecDeque, path::PathBuf};
use three_d::{
    AmbientLight, Camera, ClearState, Context, CpuMaterial, CpuMesh, DirectionalLight, FrameOutput,
    Gm, Indices, InnerSpace, Mesh, PhysicalMaterial, Positions, Srgba, Vec3, Window,
    WindowSettings, degrees, vec3,
};

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

fn classify_vertices(mesh: &mut PolygonMesh) -> VProperty<Srgba> {
    let points = mesh.points();
    // Properties must be borrowed before they can be used.
    let points = points.try_borrow().unwrap();
    // Find vertices with the most extreme coordinate values of points.
    let mut outer = mesh
        .vertices()
        .fold(None, |outer: Option<[VH; 6]>, v| match outer {
            Some(mut outer) => {
                let pos = points[v];
                if pos.x < points[outer[0]].x {
                    outer[0] = v;
                }
                if pos.x > points[outer[1]].x {
                    outer[1] = v;
                }
                if pos.y < points[outer[2]].y {
                    outer[2] = v;
                }
                if pos.y > points[outer[3]].y {
                    outer[3] = v;
                }
                if pos.z < points[outer[4]].z {
                    outer[4] = v;
                }
                if pos.z > points[outer[5]].z {
                    outer[5] = v;
                }
                Some(outer)
            }
            None => Some([v; 6]),
        })
        .unwrap()
        .to_vec();
    outer.sort();
    outer.dedup();
    let num_outer = outer.len();
    // Create vertex property to keep track of group index.
    let mut group = mesh.create_vertex_prop::<Option<usize>>(None);
    // borrow the property to be able to write to it.
    let mut group = group.try_borrow_mut().unwrap();
    // Setup the initial seeds.
    let mut vstack = VecDeque::<(usize, VH)>::with_capacity(mesh.num_vertices());
    vstack.extend(outer.into_iter().enumerate());
    // Start fast marching to classify the vertices.
    while let Some((i, v)) = vstack.pop_front() {
        if group[v].is_some() {
            // This vertex has already been visited.
            continue;
        }
        group[v] = Some(i); // Claim this vetex.
        vstack.extend(mesh.vv_ccw_iter(v).map(|v| (i, v))); // Push it's neighbors.
    }
    // Assign colors by group.
    const LEGEND: [Srgba; 6] = [
        Srgba::RED,
        Srgba::GREEN,
        Srgba::BLUE,
        Srgba::new_opaque(255, 255, 0),
        Srgba::new_opaque(255, 0, 255),
        Srgba::new_opaque(0, 255, 255),
    ];
    let legend = LEGEND.iter().take(num_outer).collect::<Vec<_>>();
    let mut colors = mesh.create_vertex_prop(Srgba::WHITE);
    {
        // Borrow colors to write into it.
        let mut colors = colors.try_borrow_mut().unwrap();
        for v in mesh.vertices() {
            let i = group[v].unwrap();
            colors[v] = *legend[i];
        }
        // The borrowed reference is dropped, allowing us to return the property from the function.
    }
    colors
}

fn mesh_view_with_colors(
    mesh: &PolygonMesh,
    context: &Context,
    colors: &VProperty<Srgba>,
) -> Gm<Mesh, PhysicalMaterial> {
    let points = mesh.points();
    let points = points.try_borrow().expect("Cannot borrow points");
    let vnormals = mesh.vertex_normals().unwrap();
    let vnormals = vnormals.try_borrow().unwrap();
    let colors = colors.try_borrow().unwrap();
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
        colors: Some(colors.to_vec()),
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

fn main() {
    // Create the scene.
    let window = Window::new(WindowSettings {
        title: "Example: Using the Property System for Colors".to_string(),
        min_size: (512, 256),
        ..Default::default()
    })
    .unwrap();
    let context = window.gl();
    // Import the mesh.
    let path = {
        let mut path = PathBuf::new();
        path.push(env!("CARGO_MANIFEST_DIR"));
        path.push("assets");
        path.push("bunny_large.obj");
        path
    };
    let mut mesh = PolygonMesh::load_obj(path).unwrap();
    mesh.update_vertex_normals_accurate().unwrap();
    let mut mesh = mesh; // Immutable.
    let (min, max) = bounds(&mesh);
    // Setup the camera and the controls and lights.
    let target = 0.5 * (min + max);
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
    // Render the mesh.
    let colors = classify_vertices(&mut mesh);
    let mview = mesh_view_with_colors(&mesh, &context, &colors);
    window.render_loop(move |mut frame_input| {
        let mut redraw = frame_input.first_frame;
        redraw |= camera.set_viewport(frame_input.viewport);
        redraw |= control.handle_events(&mut camera, &mut frame_input.events);
        if redraw {
            frame_input
                .screen()
                .clear(ClearState::color_and_depth(0.1, 0.1, 0.1, 1.0, 1.0))
                .render(&camera, &mview, &[&ambient, &directional0, &directional1]);
        }
        FrameOutput {
            swap_buffers: redraw,
            ..Default::default()
        }
    });
}
