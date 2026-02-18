//! CSG failure-scenario diagnostics.
//!
//! Reproduces current known failure modes with concrete metrics:
//! 1) Invalid operands (inward-wound cube + open sphere)
//! 2) Analytical overlapping-cube baseline (difference currently fails)
//!
//! Run with:
//! `cargo run -p cfd-mesh --example csg_failure_scenarios --features csg`

use std::f64::consts::{PI, TAU};

use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean};
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::topology::{manifold, orientation};

const EPS_VOLUME: Real = 0.1;

#[derive(Clone, Copy)]
struct MeshSummary {
    faces: usize,
    vertices: usize,
    area: Real,
    volume: Real,
    boundary_edges: usize,
    non_manifold_edges: usize,
    orientation_ok: bool,
    failing_faces: usize,
    total_faces: usize,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Failure Scenarios");
    println!("=================================================================");
    println!();

    scenario_invalid_operands()?;
    scenario_analytical_cubes()?;

    println!("=================================================================");
    println!("  Complete");
    println!("=================================================================");

    Ok(())
}

fn scenario_invalid_operands() -> Result<(), Box<dyn std::error::Error>> {
    println!("--- Scenario 1: Invalid Operands (Inward Cube + Open Sphere) ---");

    let mut pool = VertexPool::default_millifluidic();

    let mut cube = generate_cube_outward(1.0, Point3r::origin(), &mut pool, RegionId::new(1));
    for face in &mut cube {
        face.flip(); // intentionally inward
    }

    let sphere_open = generate_uv_sphere(
        Point3r::new(0.5, 0.5, 0.5),
        0.75,
        32,
        16,
        false, // intentionally open at poles
        &mut pool,
        RegionId::new(2),
    );

    let cube_summary = summarize_faces(&pool, &cube);
    let sphere_summary = summarize_faces(&pool, &sphere_open);

    print_summary("inward_cube", cube_summary);
    print_summary("open_sphere", sphere_summary);

    run_invalid_op(
        "invalid_union",
        BooleanOp::Union,
        &cube,
        &sphere_open,
        &mut pool,
        cube_summary,
        sphere_summary,
    );
    run_invalid_op(
        "invalid_intersection",
        BooleanOp::Intersection,
        &cube,
        &sphere_open,
        &mut pool,
        cube_summary,
        sphere_summary,
    );
    run_invalid_op(
        "invalid_difference",
        BooleanOp::Difference,
        &cube,
        &sphere_open,
        &mut pool,
        cube_summary,
        sphere_summary,
    );

    println!();
    Ok(())
}

fn scenario_analytical_cubes() -> Result<(), Box<dyn std::error::Error>> {
    println!("--- Scenario 2: Analytical Overlapping Cubes ---");
    println!("  Expected volumes: union=15, intersection=1, difference=7");

    let mut pool = VertexPool::default_millifluidic();
    let cube_a = generate_cube_outward(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(10));
    let cube_b = generate_cube_outward(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(20));

    run_analytic_op("union", BooleanOp::Union, 15.0, &cube_a, &cube_b, &mut pool)?;
    run_analytic_op(
        "intersection",
        BooleanOp::Intersection,
        1.0,
        &cube_a,
        &cube_b,
        &mut pool,
    )?;
    run_analytic_op(
        "difference",
        BooleanOp::Difference,
        7.0,
        &cube_a,
        &cube_b,
        &mut pool,
    )?;

    println!();
    Ok(())
}

fn run_analytic_op(
    label: &str,
    op: BooleanOp,
    expected_volume: Real,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
) -> Result<(), Box<dyn std::error::Error>> {
    let result = csg_boolean(op, faces_a, faces_b, pool)?;
    let summary = summarize_faces(pool, &result);
    let err = (summary.volume - expected_volume).abs();
    let status = if err <= EPS_VOLUME { "PASS" } else { "FAIL" };

    print_summary(label, summary);
    println!(
        "  expected_volume={:.4}, error={:.4} -> {}",
        expected_volume, err, status
    );

    Ok(())
}

fn run_invalid_op(
    label: &str,
    op: BooleanOp,
    faces_a: &[FaceData],
    faces_b: &[FaceData],
    pool: &mut VertexPool,
    a_summary: MeshSummary,
    b_summary: MeshSummary,
) {
    match csg_boolean(op, faces_a, faces_b, pool) {
        Ok(result) => {
            let result_summary = summarize_faces(pool, &result);
            print_summary(label, result_summary);
            println!(
                "  Failure markers: operand watertight={} / {} | result watertight={} | volume={:.4}",
                is_watertight(a_summary),
                is_watertight(b_summary),
                is_watertight(result_summary),
                result_summary.volume
            );
        }
        Err(err) => {
            println!("  {label} failed (expected in this scenario): {err}");
        }
    }
}

fn summarize_faces(pool: &VertexPool, faces: &[FaceData]) -> MeshSummary {
    let mut mesh = IndexedMesh::new();
    mesh.vertices = clone_pool(pool);
    for face in faces {
        mesh.faces.push(*face);
    }
    mesh.rebuild_edges();

    let area = mesh.surface_area();
    let volume = mesh.signed_volume();
    let edges = EdgeStore::from_face_store(&mesh.faces);
    let manifold_report = manifold::check_manifold(&edges);
    let orientation_ok = orientation::check_orientation(&mesh.faces, &edges).is_ok();
    let quality = mesh.quality_report();

    MeshSummary {
        faces: mesh.face_count(),
        vertices: mesh.vertex_count(),
        area,
        volume,
        boundary_edges: manifold_report.boundary_edges,
        non_manifold_edges: manifold_report.non_manifold_edges,
        orientation_ok,
        failing_faces: quality.failing_faces,
        total_faces: quality.total_faces,
    }
}

fn print_summary(label: &str, summary: MeshSummary) {
    println!(
        "  {:<14} faces={:<6} verts={:<6} area={:>9.4} vol={:>9.4} boundary={} nonmanifold={} orient_ok={} quality={}/{}",
        label,
        summary.faces,
        summary.vertices,
        summary.area,
        summary.volume,
        summary.boundary_edges,
        summary.non_manifold_edges,
        summary.orientation_ok,
        summary.failing_faces,
        summary.total_faces,
    );
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
