//! Criterion performance benchmarks for the CSG arrangement pipeline.
//!
//! ## Benchmarks
//!
//! | Benchmark | What it measures |
//! |-----------|-----------------|
//! | `gwn_linear_small` | Linear GWN for 12-face cube (baseline) |
//! | `gwn_linear_large` | Linear GWN for 10 000-face sphere (~10k faces) |
//! | `gwn_bvh_large` | BVH-accelerated GWN for 10 000-face sphere |
//! | `csg_union_cube_cube` | Full CSG union of two overlapping unit cubes |
//! | `csg_intersection_cylinders` | CSG intersection of two 64-segment cylinders |
//! | `detect_self_intersect_flat` | Self-intersection scan on a flat quad mesh |
//!
//! ## Running
//!
//! ```bash
//! cargo bench -p cfd-mesh
//! ```
//!
//! ## Expectations
//!
//! - `gwn_bvh_large` should be at least 5× faster than `gwn_linear_large` for
//!   meshes with ≥ 1000 faces (from the O(log n) vs O(n) theory).
//! - `csg_union_cube_cube` should complete in < 10 ms (regression guard).

use criterion::{black_box, criterion_group, criterion_main, Criterion};

use cfd_mesh::application::csg::arrangement::classify::{
    classify_fragment_prepared, prepare_classification_faces,
};
use cfd_mesh::application::csg::arrangement::gwn::{gwn, prepare_classification_faces as pcf};
use cfd_mesh::application::csg::arrangement::gwn_bvh::{gwn_bvh, prepare_bvh_mesh};
use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::application::csg::detect_self_intersect::detect_self_intersections;
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::{Cube, Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::infrastructure::storage::face_store::FaceData;
use cfd_mesh::infrastructure::storage::vertex_pool::VertexPool;

// ── Helper builders ────────────────────────────────────────────────────────────

fn unit_cube_faces() -> (VertexPool, Vec<FaceData>) {
    let mut pool = VertexPool::default_millifluidic();
    let n = nalgebra::Vector3::zeros();
    let s = 0.5_f64;
    let mut v = |x, y, z| pool.insert_or_weld(Point3r::new(x, y, z), n);
    let c000 = v(-s, -s, -s);
    let c100 = v(s, -s, -s);
    let c010 = v(-s, s, -s);
    let c110 = v(s, s, -s);
    let c001 = v(-s, -s, s);
    let c101 = v(s, -s, s);
    let c011 = v(-s, s, s);
    let c111 = v(s, s, s);
    let f = FaceData::untagged;
    let faces = vec![
        f(c000, c010, c110),
        f(c000, c110, c100),
        f(c001, c101, c111),
        f(c001, c111, c011),
        f(c000, c001, c011),
        f(c000, c011, c010),
        f(c100, c110, c111),
        f(c100, c111, c101),
        f(c000, c100, c101),
        f(c000, c101, c001),
        f(c010, c011, c111),
        f(c010, c111, c110),
    ];
    (pool, faces)
}

fn build_sphere_faces(rings: u32, sectors: u32) -> (VertexPool, Vec<FaceData>) {
    use cfd_mesh::domain::mesh::IndexedMesh;
    let mesh = UvSphere { center: Point3r::origin(), radius: 1.0, rings, sectors }
        .build()
        .expect("sphere build");
    let mut pool = VertexPool::default_millifluidic();
    let mut id_map: std::collections::HashMap<
        cfd_mesh::domain::core::index::VertexId,
        cfd_mesh::domain::core::index::VertexId,
    > = std::collections::HashMap::new();
    for (old, _) in mesh.vertices.iter() {
        let pos = *mesh.vertices.position(old);
        let nrm = *mesh.vertices.normal(old);
        id_map.insert(old, pool.insert_or_weld(pos, nrm));
    }
    let faces = mesh
        .faces
        .iter()
        .map(|f| FaceData { vertices: f.vertices.map(|v| id_map[&v]), region: f.region })
        .collect();
    (pool, faces)
}

// ── GWN benchmarks ─────────────────────────────────────────────────────────────

fn bench_gwn_linear_small(c: &mut Criterion) {
    let (pool, faces) = unit_cube_faces();
    let q = Point3r::new(0.0, 0.0, 0.0);
    c.bench_function("gwn_linear_small_12f", |b| {
        b.iter(|| gwn::<f64>(black_box(&q), black_box(&faces), black_box(&pool)))
    });
}

fn bench_gwn_linear_large(c: &mut Criterion) {
    // ~2400 faces for rings=40, sectors=60
    let (pool, faces) = build_sphere_faces(40, 60);
    let q = Point3r::new(0.0, 0.0, 0.0); // interior query
    c.bench_function("gwn_linear_2400f", |b| {
        b.iter(|| gwn::<f64>(black_box(&q), black_box(&faces), black_box(&pool)))
    });
}

fn bench_gwn_bvh_large(c: &mut Criterion) {
    let (pool, faces) = build_sphere_faces(40, 60);
    let prepared = prepare_bvh_mesh(&faces, &pool);
    let q = Point3r::new(0.0, 0.0, 0.0);
    c.bench_function("gwn_bvh_2400f", |b| {
        b.iter(|| gwn_bvh(black_box(&q), black_box(&prepared), black_box(0.01)))
    });
}

// ── Prepared classification benchmarks ────────────────────────────────────────

fn bench_classify_prepared(c: &mut Criterion) {
    let (pool, faces) = build_sphere_faces(40, 60);
    let prepared = pcf(&faces, &pool);
    let q = Point3r::new(0.0, 0.0, 0.0);
    let n = nalgebra::Vector3::new(0.0, 0.0, 1.0);
    c.bench_function("classify_prepared_2400f", |b| {
        b.iter(|| classify_fragment_prepared(black_box(&q), black_box(&n), black_box(&prepared)))
    });
}

// ── Full CSG benchmarks ────────────────────────────────────────────────────────

fn bench_csg_union_cube_cube(c: &mut Criterion) {
    let cube_a = Cube { origin: Point3r::new(-1.0, -1.0, -1.0), width: 2.0, height: 2.0, depth: 2.0 }
        .build().unwrap();
    let cube_b = Cube { origin: Point3r::new(-0.5, -0.5, -0.5), width: 2.0, height: 2.0, depth: 2.0 }
        .build().unwrap();
    c.bench_function("csg_union_cube_cube", |b| {
        b.iter(|| {
            csg_boolean_indexed(BooleanOp::Union, black_box(&cube_a), black_box(&cube_b)).ok()
        })
    });
}

fn bench_csg_intersection_cylinders(c: &mut Criterion) {
    let cyl_a = Cylinder {
        base_center: Point3r::new(-2.0, 0.0, 0.0),
        radius: 1.0,
        height: 4.0,
        segments: 32,
    }.build().unwrap();
    let cyl_b = Cylinder {
        base_center: Point3r::new(0.0, -2.0, 0.0),
        radius: 1.0,
        height: 4.0,
        segments: 32,
    }.build().unwrap();
    c.bench_function("csg_intersection_cylinders_32seg", |b| {
        b.iter(|| {
            csg_boolean_indexed(BooleanOp::Intersection, black_box(&cyl_a), black_box(&cyl_b)).ok()
        })
    });
}

// ── Self-intersection scan benchmark ──────────────────────────────────────────

fn bench_detect_self_intersect_flat(c: &mut Criterion) {
    // Build a flat 10×10 quad grid of triangles (200 triangles) — no self-intersections.
    let mut pool = VertexPool::default_millifluidic();
    let n = nalgebra::Vector3::zeros();
    let mut faces = Vec::new();
    let n_side = 10_usize;
    let mut verts = vec![vec![cfd_mesh::domain::core::index::VertexId::new(0); n_side + 1]; n_side + 1];
    for i in 0..=n_side {
        for j in 0..=n_side {
            verts[i][j] = pool.insert_or_weld(
                Point3r::new(i as f64, j as f64, 0.0),
                n,
            );
        }
    }
    for i in 0..n_side {
        for j in 0..n_side {
            faces.push(FaceData::untagged(verts[i][j], verts[i+1][j], verts[i+1][j+1]));
            faces.push(FaceData::untagged(verts[i][j], verts[i+1][j+1], verts[i][j+1]));
        }
    }
    c.bench_function("detect_self_intersect_200tri_flat", |b| {
        b.iter(|| detect_self_intersections(black_box(&faces), black_box(&pool)))
    });
}

// ── Criterion groups ───────────────────────────────────────────────────────────

criterion_group!(
    gwn_benches,
    bench_gwn_linear_small,
    bench_gwn_linear_large,
    bench_gwn_bvh_large,
    bench_classify_prepared,
);

criterion_group!(
    csg_benches,
    bench_csg_union_cube_cube,
    bench_csg_intersection_cylinders,
);

criterion_group!(detect_benches, bench_detect_self_intersect_flat);

criterion_main!(gwn_benches, csg_benches, detect_benches);
