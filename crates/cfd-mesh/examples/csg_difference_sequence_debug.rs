//! Difference-sequence debugger for CSG BSP operations.
//!
//! Compares:
//! - Current library `BooleanOp::Difference`
//! - Manual sequence with an additional `b.invert()` before merge
//!
//! Uses overlapping cubes with analytical expected volume 7.0 mm³.
//!
//! Run with:
//! `cargo run -p cfd-mesh --example csg_difference_sequence_debug --features csg`

use cfd_mesh::IndexedMesh;
use cfd_mesh::core::index::RegionId;
use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean};
use cfd_mesh::csg::bsp::BspNode;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::storage::face_store::FaceData;
use cfd_mesh::storage::vertex_pool::VertexPool;
use cfd_mesh::topology::{manifold, orientation};

const EXPECTED_DIFFERENCE_VOLUME: Real = 7.0;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Difference Sequence Debugger");
    println!("=================================================================");
    println!("  Geometry: cube A [0,2]^3, cube B [1,3]^3");
    println!("  Expected A\\B volume: {:.1} mm³", EXPECTED_DIFFERENCE_VOLUME);
    println!();

    let current = run_current_difference()?;
    let patched = run_manual_difference_with_extra_invert()?;

    println!("--- Comparison ---");
    print_result("current_boolean_difference", current);
    print_result("manual_with_extra_invert", patched);

    let current_err = (current.volume - EXPECTED_DIFFERENCE_VOLUME).abs();
    let patched_err = (patched.volume - EXPECTED_DIFFERENCE_VOLUME).abs();

    println!();
    println!("  error_current={:.4}, error_manual={:.4}", current_err, patched_err);
    if patched_err < current_err {
        println!("  Signal: manual sequence is closer to analytical volume.");
    } else {
        println!("  Signal: no improvement from manual sequence on this build.");
    }

    Ok(())
}

#[derive(Clone, Copy)]
struct ResultSummary {
    faces: usize,
    vertices: usize,
    area: Real,
    volume: Real,
    boundary_edges: usize,
    non_manifold_edges: usize,
    orientation_ok: bool,
}

fn run_current_difference() -> Result<ResultSummary, Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let a = generate_cube(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(1));
    let b = generate_cube(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(2));

    let faces = csg_boolean(BooleanOp::Difference, &a, &b, &mut pool)?;
    Ok(summarize(&pool, &faces))
}

fn run_manual_difference_with_extra_invert() -> Result<ResultSummary, Box<dyn std::error::Error>> {
    let mut pool = VertexPool::default_millifluidic();
    let a_faces = generate_cube(2.0, Point3r::new(0.0, 0.0, 0.0), &mut pool, RegionId::new(1));
    let b_faces = generate_cube(2.0, Point3r::new(1.0, 1.0, 1.0), &mut pool, RegionId::new(2));

    let mut a = BspNode::build(&a_faces, &mut pool);
    let mut b = BspNode::build(&b_faces, &mut pool);

    a.invert();
    a.clip_to(&b, &mut pool);
    b.clip_to(&a, &mut pool);
    b.invert();
    b.clip_to(&a, &mut pool);
    b.invert(); // extra invert for debugging comparison

    let b_remaining = b.all_faces();
    a.add_faces(&b_remaining, &mut pool);
    a.invert();

    let result = a.all_faces();
    Ok(summarize(&pool, &result))
}

fn summarize(pool: &VertexPool, faces: &[FaceData]) -> ResultSummary {
    let mut mesh = IndexedMesh::new();
    mesh.vertices = clone_pool(pool);
    for face in faces {
        mesh.faces.push(*face);
    }
    mesh.rebuild_edges();

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let manifold_report = manifold::check_manifold(&edges);
    let orientation_ok = orientation::check_orientation(&mesh.faces, &edges).is_ok();

    ResultSummary {
        faces: mesh.face_count(),
        vertices: mesh.vertex_count(),
        area: mesh.surface_area(),
        volume: mesh.signed_volume(),
        boundary_edges: manifold_report.boundary_edges,
        non_manifold_edges: manifold_report.non_manifold_edges,
        orientation_ok,
    }
}

fn print_result(label: &str, summary: ResultSummary) {
    println!(
        "  {:<28} faces={:<6} verts={:<6} area={:>9.4} vol={:>9.4} boundary={} nonmanifold={} orient_ok={}",
        label,
        summary.faces,
        summary.vertices,
        summary.area,
        summary.volume,
        summary.boundary_edges,
        summary.non_manifold_edges,
        summary.orientation_ok,
    );
}

fn clone_pool(pool: &VertexPool) -> VertexPool {
    let mut cloned = VertexPool::default_millifluidic();
    for (_, v) in pool.iter() {
        cloned.insert_unique(v.position, v.normal);
    }
    cloned
}

fn generate_cube(
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
