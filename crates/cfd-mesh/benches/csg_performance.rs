//! Criterion performance benchmarks for the CSG arrangement pipeline.
use criterion::{black_box, criterion_group, criterion_main, Criterion};

use cfd_mesh::application::csg::arrangement::classify::{
    classify_fragment_prepared, gwn, gwn_bvh, prepare_bvh_mesh, prepare_classification_faces,
};
use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::application::csg::detect_self_intersect::detect_self_intersections;
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::{Cube, Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::infrastructure::storage::face_store::FaceData;
use cfd_mesh::infrastructure::storage::vertex_pool::VertexPool;

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
    (
        pool,
        vec![
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
        ],
    )
}

fn build_sphere_faces(stacks: usize, segments: usize) -> (VertexPool, Vec<FaceData>) {
    let mesh = UvSphere {
        center: Point3r::origin(),
        radius: 1.0,
        stacks,
        segments,
    }
    .build()
    .expect("sphere");
    let mut pool = VertexPool::default_millifluidic();
    let mut id_map = std::collections::HashMap::new();
    for (old, _) in mesh.vertices.iter() {
        id_map.insert(
            old,
            pool.insert_or_weld(*mesh.vertices.position(old), *mesh.vertices.normal(old)),
        );
    }
    let faces = mesh
        .faces
        .iter()
        .map(|f| FaceData {
            vertices: f.vertices.map(|v| id_map[&v]),
            region: f.region,
        })
        .collect();
    (pool, faces)
}

fn bench_gwn_linear_small(c: &mut Criterion) {
    let (pool, faces) = unit_cube_faces();
    let q = Point3r::new(0.0, 0.0, 0.0);
    c.bench_function("gwn_linear_12f", |b| {
        b.iter(|| gwn::<f64>(black_box(&q), black_box(&faces), black_box(&pool)))
    });
}

fn bench_gwn_linear_large(c: &mut Criterion) {
    let (pool, faces) = build_sphere_faces(40, 60);
    let q = Point3r::new(0.0, 0.0, 0.0);
    c.bench_function("gwn_linear_2400f", |b| {
        b.iter(|| gwn::<f64>(black_box(&q), black_box(&faces), black_box(&pool)))
    });
}

fn bench_gwn_bvh_large(c: &mut Criterion) {
    let (pool, faces) = build_sphere_faces(40, 60);
    // prepare_bvh_mesh takes &[PreparedFace] (pre-computed via prepare_classification_faces).
    let prep_faces = prepare_classification_faces(&faces, &pool);
    let bvh = prepare_bvh_mesh(&prep_faces).expect("non-empty sphere should build BVH");
    let q = Point3r::new(0.0, 0.0, 0.0);
    c.bench_function("gwn_bvh_2400f", |b| {
        b.iter(|| gwn_bvh(black_box(&q), black_box(&bvh), black_box(0.01)))
    });
}

fn bench_classify_prepared(c: &mut Criterion) {
    let (pool, faces) = build_sphere_faces(40, 60);
    let prepared = prepare_classification_faces(&faces, &pool);
    let q = Point3r::new(0.0, 0.0, 0.0);
    let n = nalgebra::Vector3::new(0.0_f64, 0.0, 1.0);
    c.bench_function("classify_prepared_2400f", |b| {
        b.iter(|| classify_fragment_prepared(black_box(&q), black_box(&n), black_box(&prepared)))
    });
}

fn bench_csg_union_cube_cube(c: &mut Criterion) {
    let cube_a = Cube {
        origin: Point3r::new(-1.0, -1.0, -1.0),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();
    let cube_b = Cube {
        origin: Point3r::new(-0.5, -0.5, -0.5),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()
    .unwrap();
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
    }
    .build()
    .unwrap();
    let cyl_b = Cylinder {
        base_center: Point3r::new(0.0, -2.0, 0.0),
        radius: 1.0,
        height: 4.0,
        segments: 32,
    }
    .build()
    .unwrap();
    c.bench_function("csg_intersection_cyl32", |b| {
        b.iter(|| {
            csg_boolean_indexed(
                BooleanOp::Intersection,
                black_box(&cyl_a),
                black_box(&cyl_b),
            )
            .ok()
        })
    });
}

fn bench_detect_self_intersect_flat(c: &mut Criterion) {
    let mut pool = VertexPool::default_millifluidic();
    let n = nalgebra::Vector3::zeros();
    let ns = 10_usize;
    let mut verts = vec![vec![cfd_mesh::domain::core::index::VertexId::new(0); ns + 1]; ns + 1];
    for i in 0..=ns {
        for j in 0..=ns {
            verts[i][j] = pool.insert_or_weld(Point3r::new(i as f64, j as f64, 0.0), n);
        }
    }
    let mut faces = Vec::new();
    for i in 0..ns {
        for j in 0..ns {
            faces.push(FaceData::untagged(
                verts[i][j],
                verts[i + 1][j],
                verts[i + 1][j + 1],
            ));
            faces.push(FaceData::untagged(
                verts[i][j],
                verts[i + 1][j + 1],
                verts[i][j + 1],
            ));
        }
    }
    c.bench_function("detect_self_intersect_200tri", |b| {
        b.iter(|| detect_self_intersections(black_box(&faces), black_box(&pool)))
    });
}

criterion_group!(
    gwn_benches,
    bench_gwn_linear_small,
    bench_gwn_linear_large,
    bench_gwn_bvh_large,
    bench_classify_prepared
);
criterion_group!(
    csg_benches,
    bench_csg_union_cube_cube,
    bench_csg_intersection_cylinders
);
criterion_group!(detect_benches, bench_detect_self_intersect_flat);
criterion_main!(gwn_benches, csg_benches, detect_benches);
