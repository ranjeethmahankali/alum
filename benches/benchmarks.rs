use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;

// Define mesh adaptor similar to examples
use alum::{HasTopology, use_glam::PolyMeshF32};

type PolygonMesh = PolyMeshF32;

// Helper function to get asset path
fn get_asset_path(filename: &str) -> PathBuf {
    let mut path = PathBuf::new();
    path.push(env!("CARGO_MANIFEST_DIR"));
    path.push("assets");
    path.push(filename);
    path
}

// OBJ Loading Benchmarks
fn bench_obj_loading(c: &mut Criterion) {
    let mut group = c.benchmark_group("obj_loading");

    // Benchmark loading the bunny mesh
    group.bench_function("load_bunny", |b| {
        let path = get_asset_path("bunny.obj");
        b.iter(|| {
            let mesh = PolygonMesh::load_obj(black_box(&path)).unwrap();
            black_box(mesh);
        });
    });

    // Benchmark loading the large bunny mesh
    group.bench_function("load_bunny_large", |b| {
        let path = get_asset_path("bunny_large.obj");
        b.iter(|| {
            let mesh = PolygonMesh::load_obj(black_box(&path)).unwrap();
            black_box(mesh);
        });
    });

    // Benchmark loading the quad bunny mesh
    group.bench_function("load_bunny_quad", |b| {
        let path = get_asset_path("bunny_quad.obj");
        b.iter(|| {
            let mesh = PolygonMesh::load_obj(black_box(&path)).unwrap();
            black_box(mesh);
        });
    });

    group.finish();
}

// Primitive Creation Benchmarks
fn bench_primitive_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("primitive_creation");

    // Benchmark creating multiple tetrahedra
    group.bench_function("tetrahedron_batch", |b| {
        b.iter(|| {
            for i in 0..100 {
                let mesh = PolygonMesh::tetrahedron(black_box(1.0 + i as f32 * 0.01)).unwrap();
                black_box(mesh);
            }
        });
    });

    // Benchmark creating multiple hexahedra (cubes)
    group.bench_function("hexahedron_batch", |b| {
        b.iter(|| {
            for i in 0..100 {
                let mesh = PolygonMesh::hexahedron(black_box(1.0 + i as f32 * 0.01)).unwrap();
                black_box(mesh);
            }
        });
    });

    // Benchmark creating multiple octahedra
    group.bench_function("octahedron_batch", |b| {
        b.iter(|| {
            for i in 0..100 {
                let mesh = PolygonMesh::octahedron(black_box(1.0 + i as f32 * 0.01)).unwrap();
                black_box(mesh);
            }
        });
    });

    // Benchmark creating multiple dodecahedra
    group.bench_function("dodecahedron_batch", |b| {
        b.iter(|| {
            for i in 0..50 {
                // Fewer iterations as dodecahedron is more complex
                let mesh = PolygonMesh::dodecahedron(black_box(1.0 + i as f32 * 0.01)).unwrap();
                black_box(mesh);
            }
        });
    });

    // Benchmark creating multiple icosahedra
    group.bench_function("icosahedron_batch", |b| {
        b.iter(|| {
            for i in 0..50 {
                // Fewer iterations as icosahedron is more complex
                let mesh = PolygonMesh::icosahedron(black_box(1.0 + i as f32 * 0.01)).unwrap();
                black_box(mesh);
            }
        });
    });

    // Benchmark creating unit boxes
    group.bench_function("unit_box_batch", |b| {
        b.iter(|| {
            for _i in 0..100 {
                let mesh = PolygonMesh::unit_box().unwrap();
                black_box(mesh);
            }
        });
    });

    group.finish();
}

// Subdivision Benchmarks
fn bench_subdivision(c: &mut Criterion) {
    let mut group = c.benchmark_group("subdivision");

    // Load bunny mesh once for all subdivision benchmarks
    let bunny_path = get_asset_path("bunny.obj");
    let base_bunny = PolygonMesh::load_obj(&bunny_path).unwrap();

    // Catmull-Clark subdivision
    group.bench_function("catmull_clark_bunny", |b| {
        b.iter(|| {
            let mut mesh = base_bunny.try_clone().unwrap();
            mesh.subdivide_catmull_clark(black_box(1), true).unwrap();
            black_box(mesh);
        });
    });

    // Loop subdivision
    group.bench_function("loop_subdivision_bunny", |b| {
        b.iter(|| {
            let mut mesh = base_bunny.try_clone().unwrap();
            mesh.subdivide_loop(black_box(1), true).unwrap();
            black_box(mesh);
        });
    });

    // Sqrt3 subdivision
    group.bench_function("sqrt3_subdivision_bunny", |b| {
        b.iter(|| {
            let mut mesh = base_bunny.try_clone().unwrap();
            mesh.subdivide_sqrt3(black_box(1), true).unwrap();
            black_box(mesh);
        });
    });

    // Multiple iterations of Catmull-Clark
    group.bench_function("catmull_clark_bunny_2_iterations", |b| {
        b.iter(|| {
            let mut mesh = base_bunny.try_clone().unwrap();
            mesh.subdivide_catmull_clark(black_box(2), true).unwrap();
            black_box(mesh);
        });
    });

    // Multiple iterations of Loop subdivision
    group.bench_function("loop_subdivision_bunny_2_iterations", |b| {
        b.iter(|| {
            let mut mesh = base_bunny.try_clone().unwrap();
            mesh.subdivide_loop(black_box(2), true).unwrap();
            black_box(mesh);
        });
    });

    // Test subdivision on simpler mesh for comparison
    let cube = PolygonMesh::unit_box().unwrap();

    group.bench_function("catmull_clark_cube", |b| {
        b.iter(|| {
            let mut mesh = cube.try_clone().unwrap();
            mesh.subdivide_catmull_clark(black_box(1), true).unwrap();
            black_box(mesh);
        });
    });

    group.bench_function("loop_subdivision_cube", |b| {
        b.iter(|| {
            let mut mesh = cube.try_clone().unwrap();
            mesh.subdivide_loop(black_box(1), true).unwrap();
            black_box(mesh);
        });
    });

    group.finish();
}

// Additional operations benchmarks
fn bench_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("operations");

    let bunny_path = get_asset_path("bunny.obj");
    let bunny = PolygonMesh::load_obj(&bunny_path).unwrap();

    // Benchmark normal calculations
    group.bench_function("update_vertex_normals_accurate", |b| {
        b.iter(|| {
            let mut mesh = bunny.try_clone().unwrap();
            mesh.update_face_normals().unwrap();
            mesh.update_vertex_normals_accurate().unwrap();
            black_box(mesh);
        });
    });

    group.bench_function("update_vertex_normals_fast", |b| {
        b.iter(|| {
            let mut mesh = bunny.try_clone().unwrap();
            mesh.update_face_normals().unwrap();
            mesh.update_vertex_normals_fast().unwrap();
            black_box(mesh);
        });
    });

    // Benchmark cloning
    group.bench_function("clone", |b| {
        b.iter(|| {
            let mesh_clone = bunny.try_clone().unwrap();
            black_box(mesh_clone);
        });
    });

    // Benchmark garbage collection
    group.bench_function("garbage_collection", |b| {
        b.iter(|| {
            let mut mesh = bunny.try_clone().unwrap();
            mesh.garbage_collection().unwrap();
            black_box(mesh);
        });
    });

    // Benchmark topology iteration
    group.bench_function("vertex_iteration", |b| {
        b.iter(|| {
            let count = bunny.vertices().count();
            black_box(count);
        });
    });

    group.bench_function("face_iteration", |b| {
        b.iter(|| {
            let count = bunny.faces().count();
            black_box(count);
        });
    });

    group.bench_function("edge_iteration", |b| {
        b.iter(|| {
            let count = bunny.edges().count();
            black_box(count);
        });
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_obj_loading,
    bench_primitive_creation,
    bench_subdivision,
    bench_operations
);
criterion_main!(benches);
