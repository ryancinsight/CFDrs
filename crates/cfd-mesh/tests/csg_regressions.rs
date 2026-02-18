#![cfg(feature = "csg")]

use std::f64::consts::{PI, TAU};

use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean};
use cfd_mesh::csg::bsp::BspNode;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::topology::{manifold, orientation};

const EPS_VOLUME: Real = 0.1;

#[derive(Clone, Copy)]
struct MeshSummary {
    volume: Real,
    boundary_edges: usize,
    non_manifold_edges: usize,
    orientation_ok: bool,
}

#[test]
fn overlapping_cubes_match_analytical_volumes() -> Result<(), Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let cube_a = generate_cube_outward(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(1));
    let cube_b = generate_cube_outward(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(2));

    let union_faces = csg_boolean(BooleanOp::Union, &cube_a, &cube_b, &mut pool)?;
    let intersection_faces = csg_boolean(BooleanOp::Intersection, &cube_a, &cube_b, &mut pool)?;
    let difference_faces = csg_boolean(BooleanOp::Difference, &cube_a, &cube_b, &mut pool)?;

    let union = summarize_faces(&pool, &union_faces);
    let intersection = summarize_faces(&pool, &intersection_faces);
    let difference = summarize_faces(&pool, &difference_faces);

    assert!((union.volume - 15.0).abs() <= EPS_VOLUME, "union volume={}", union.volume);
    assert!(
        (intersection.volume - 1.0).abs() <= EPS_VOLUME,
        "intersection volume={}",
        intersection.volume
    );
    assert!(
        (difference.volume - 7.0).abs() <= EPS_VOLUME,
        "difference volume={}",
        difference.volume
    );

    Ok(())
}

#[test]
fn invalid_operands_reproduce_failure_markers() {
    let mut pool = VertexPool::default_millifluidic();

    let mut inward_cube = generate_cube_outward(1.0, Point3r::origin(), &mut pool, RegionId::new(10));
    for face in &mut inward_cube {
        face.flip();
    }

    let open_sphere =
        generate_uv_sphere(Point3r::new(0.5, 0.5, 0.5), 0.75, 32, 16, false, &mut pool, RegionId::new(20));

    for op in [BooleanOp::Union, BooleanOp::Intersection, BooleanOp::Difference] {
        match csg_boolean(op, &inward_cube, &open_sphere, &mut pool) {
            Ok(result) => {
                let summary = summarize_faces(&pool, &result);
                let has_failure_marker = !is_watertight(summary) || summary.volume <= 0.0;
                assert!(
                    has_failure_marker,
                    "expected failure markers for {:?}, got volume={} boundary={} non_manifold={} orientation_ok={}",
                    op,
                    summary.volume,
                    summary.boundary_edges,
                    summary.non_manifold_edges,
                    summary.orientation_ok
                );
            }
            Err(_) => {}
        }
    }
}

#[test]
fn difference_sequence_manual_variant_matches_current_on_cube_case(
) -> Result<(), Box<dyn std::error::Error>> {
    let current = run_current_difference()?;
    let manual = run_manual_difference_with_extra_invert()?;

    assert!(
        (current.volume - manual.volume).abs() <= EPS_VOLUME,
        "current={} manual={}",
        current.volume,
        manual.volume
    );

    Ok(())
}

fn run_current_difference() -> Result<MeshSummary, Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let a = generate_cube_outward(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(1));
    let b = generate_cube_outward(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(2));
    let faces = csg_boolean(BooleanOp::Difference, &a, &b, &mut pool)?;
    Ok(summarize_faces(&pool, &faces))
}

fn run_manual_difference_with_extra_invert() -> Result<MeshSummary, Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let a_faces = generate_cube_outward(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(1));
    let b_faces = generate_cube_outward(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(2));

    let mut a = BspNode::build(&a_faces, &mut pool);
    let mut b = BspNode::build(&b_faces, &mut pool);

    a.invert();
    a.clip_to(&b, &mut pool);
    b.clip_to(&a, &mut pool);
    b.invert();
    b.clip_to(&a, &mut pool);
    b.invert();

    let b_remaining = b.all_faces();
    a.add_faces(&b_remaining, &mut pool);
    a.invert();

    Ok(summarize_faces(&pool, &a.all_faces()))
}

fn summarize_faces(pool: &VertexPool, faces: &[FaceData]) -> MeshSummary {
    let mut mesh = IndexedMesh::new();
    mesh.vertices = clone_pool(pool);
    for face in faces {
        mesh.faces.push(*face);
    }
    mesh.rebuild_edges();

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let manifold_report = manifold::check_manifold(&edges);
    let orientation_ok = orientation::check_orientation(&mesh.faces, &edges).is_ok();

    MeshSummary {
        volume: mesh.signed_volume(),
        boundary_edges: manifold_report.boundary_edges,
        non_manifold_edges: manifold_report.non_manifold_edges,
        orientation_ok,
    }
}

fn is_watertight(summary: MeshSummary) -> bool {
    summary.boundary_edges == 0 && summary.non_manifold_edges == 0 && summary.orientation_ok
}

fn clone_pool(pool: &VertexPool) -> VertexPool {
    let mut cloned = VertexPool::default_millifluidic();
    for (_, v) in pool.iter() {
        cloned.insert_unique(v.position, v.normal);
    }
    cloned
}

fn generate_cube_outward(
    size: Real,
    origin: Point3r,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    let s = size;
    let o = origin;
    let mut faces = Vec::with_capacity(12);

    let p000 = Point3r::new(o.x, o.y, o.z);
    let p100 = Point3r::new(o.x + s, o.y, o.z);
    let p110 = Point3r::new(o.x + s, o.y + s, o.z);
    let p010 = Point3r::new(o.x, o.y + s, o.z);
    let p001 = Point3r::new(o.x, o.y, o.z + s);
    let p101 = Point3r::new(o.x + s, o.y, o.z + s);
    let p111 = Point3r::new(o.x + s, o.y + s, o.z + s);
    let p011 = Point3r::new(o.x, o.y + s, o.z + s);

    let mut add_quad = |p0: Point3r, p1: Point3r, p2: Point3r, p3: Point3r, normal: Vector3r| {
        let v0 = pool.insert_or_weld(p0, normal);
        let v1 = pool.insert_or_weld(p1, normal);
        let v2 = pool.insert_or_weld(p2, normal);
        let v3 = pool.insert_or_weld(p3, normal);
        faces.push(FaceData::new(v0, v1, v2, region));
        faces.push(FaceData::new(v0, v2, v3, region));
    };

    add_quad(p000, p010, p110, p100, -Vector3r::z());
    add_quad(p001, p101, p111, p011, Vector3r::z());
    add_quad(p000, p100, p101, p001, -Vector3r::y());
    add_quad(p010, p011, p111, p110, Vector3r::y());
    add_quad(p000, p001, p011, p010, -Vector3r::x());
    add_quad(p100, p110, p111, p101, Vector3r::x());

    faces
}

fn generate_uv_sphere(
    center: Point3r,
    radius: Real,
    segments: usize,
    stacks: usize,
    include_caps: bool,
    pool: &mut VertexPool,
    region: RegionId,
) -> Vec<FaceData> {
    let mut faces = Vec::with_capacity(segments * stacks * 2);

    let vertex_at = |theta: Real, phi: Real| -> (Point3r, Vector3r) {
        let sin_phi = phi.sin();
        let cos_phi = phi.cos();
        let sin_theta = theta.sin();
        let cos_theta = theta.cos();

        let normal = Vector3r::new(sin_phi * cos_theta, cos_phi, sin_phi * sin_theta);
        let position = center + normal * radius;
        (position, normal)
    };

    for i in 0..segments {
        for j in 0..stacks {
            let t0 = i as Real / segments as Real;
            let t1 = (i + 1) as Real / segments as Real;
            let p0 = j as Real / stacks as Real;
            let p1 = (j + 1) as Real / stacks as Real;

            let theta0 = t0 * TAU;
            let theta1 = t1 * TAU;
            let phi0 = p0 * PI;
            let phi1 = p1 * PI;

            let (pos00, n00) = vertex_at(theta0, phi0);
            let (pos10, n10) = vertex_at(theta1, phi0);
            let (pos11, n11) = vertex_at(theta1, phi1);
            let (pos01, n01) = vertex_at(theta0, phi1);

            let v00 = pool.insert_or_weld(pos00, n00);
            let v10 = pool.insert_or_weld(pos10, n10);
            let v11 = pool.insert_or_weld(pos11, n11);
            let v01 = pool.insert_or_weld(pos01, n01);

            if j == 0 {
                if include_caps {
                    faces.push(FaceData::new(v10, v11, v01, region));
                }
            } else if j == stacks - 1 {
                if include_caps {
                    faces.push(FaceData::new(v00, v10, v01, region));
                }
            } else {
                faces.push(FaceData::new(v00, v10, v11, region));
                faces.push(FaceData::new(v00, v11, v01, region));
            }
        }
    }

    faces
}