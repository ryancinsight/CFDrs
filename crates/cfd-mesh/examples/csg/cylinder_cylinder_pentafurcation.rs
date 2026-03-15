//! CSG N-Way Cylinder Pentafurcation
//!
//! One **trunk** cylinder along +X ends at a junction from which five
//! **branches** diverge in the XY plane radially:
//!
//! - Branch 1: at +60° from the +X axis
//! - Branch 2: at +30° from the +X axis
//! - Branch 3: at +15° from the +X axis (near-forward)
//! - Branch 4: at -30° from the +X axis
//! - Branch 5: at -60° from the +X axis
//!
//! This example exercises the canonical indexed N-way Boolean path on a dense
//! five-branch junction.  The center branch uses 15° (not 0°) to avoid
//! coincident surfaces with the trunk lateral barrel.
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_pentafurcation
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

use cfd_mesh::application::csg::boolean::BooleanOp;
use cfd_mesh::application::csg::CsgNode;
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::domain::topology::AdjacencyGraph;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, IndexedMesh};

// ── Geometry parameters ────────────────────────────────────────────────────────

const R: f64 = 0.5;
const H_TRUNK: f64 = 3.0;
const H_BRANCH: f64 = 3.0;
const EPS: f64 = R * 0.10;
const SEGS: usize = 64;

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG N-Way (Pentafurcation)");
    println!("  Outer branches ±60°, ±30° | Near-forward 15° branch");
    println!("  Unified N-Way Exact Union");
    println!("  r={R} mm  h_trunk={H_TRUNK} mm  h_branch={H_BRANCH} mm");
    println!("=================================================================");
    println!();

    let v_a = std::f64::consts::PI * R * R * (H_TRUNK + EPS);
    let v_b = std::f64::consts::PI * R * R * H_BRANCH;
    let v_naive = v_a + 5.0 * v_b;

    println!("  V_A (trunk+ε)       = {v_a:.4} mm³");
    println!("  V_B (x 5 branches)  = {v_b:.4} mm³ (each)");
    println!("  V_naive (sum)       = {v_naive:.4} mm³");
    println!("  Actual union        < V_naive");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let meshes = build_pentafurcation()?;
    println!("  Meshes generated  ({} ms)", t_build.elapsed().as_millis());

    for op in [
        BooleanOp::Union,
        BooleanOp::Intersection,
        BooleanOp::Difference,
    ] {
        let op_name = match op {
            BooleanOp::Union => "Union",
            BooleanOp::Intersection => "Intersection",
            BooleanOp::Difference => "Difference",
        };

        let t0 = Instant::now();
        let mut accumulated = cfd_mesh::application::csg::boolean::indexed::csg_boolean_nary(op, &meshes)?;
        let ms = t0.elapsed().as_millis();

        report(
            &format!("N-Way Pentafurcation {}", op_name),
            &mut accumulated,
            v_naive,
            ms,
            2,
        );

        let stl_name = format!(
            "cylinder_cylinder_pentafurcation_{}.stl",
            op_name.to_lowercase()
        );
        if accumulated.face_count() > 0 {
            write_stl(&accumulated, &out_dir.join(&stl_name))?;
            println!("  STL: outputs/csg/{stl_name}");
        } else {
            println!("  STL skipped: 0 faces generated for {}", op_name);
        }
        println!("-----------------------------------------------------------------");
    }
    println!("=================================================================");

    Ok(())
}

fn build_pentafurcation() -> Result<Vec<IndexedMesh>, Box<dyn std::error::Error>> {
    let mut out = Vec::new();

    // Trunk
    let raw = Cylinder {
        base_center: Point3r::new(0.0, 0.0, 0.0),
        radius: R,
        height: H_TRUNK + EPS,
        segments: SEGS,
    }
    .build()?;
    let rot =
        UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -std::f64::consts::FRAC_PI_2);
    let iso = Isometry3::from_parts(Translation3::new(-H_TRUNK, 0.0, 0.0), rot);
    out.push(
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(Box::new(raw))),
            iso,
        }
        .evaluate()?,
    );

    // Branches at +60°, +30°, +15°, -30°, -60°
    for &angle_deg in &[60.0, 30.0, 15.0, -30.0, -60.0_f64] {
        out.push(make_branch_planar(angle_deg.to_radians())?);
    }

    Ok(out)
}

fn make_branch_planar(angle_from_x: f64) -> Result<IndexedMesh, Box<dyn std::error::Error>> {
    let raw = Cylinder {
        base_center: Point3r::new(0.0, 0.0, 0.0),
        radius: R,
        height: H_BRANCH,
        segments: SEGS,
    }
    .build()?;
    let rot = UnitQuaternion::<Real>::from_axis_angle(
        &Vector3::z_axis(),
        angle_from_x - std::f64::consts::FRAC_PI_2,
    );
    let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
    Ok(CsgNode::Transform {
        node: Box::new(CsgNode::Leaf(Box::new(raw))),
        iso,
    }
    .evaluate()?)
}

// ── Report helpers ─────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, v_naive: f64, ms: u128, expected_chi: i64) {
    let vol = mesh.signed_volume();
    let n = analyze_normals(mesh);
    let positive = vol > 0.0;
    let below_naive = vol < v_naive;

    mesh.rebuild_edges();
    let wt = check_watertight(&mesh.vertices, &mesh.faces, mesh.edges_ref().unwrap());
    let adj = AdjacencyGraph::build(&mesh.faces, mesh.edges_ref().unwrap());
    let n_comps = connected_components(&mesh.faces, &adj).len();

    let chi_ok = wt.euler_characteristic == Some(expected_chi);
    let comps_ok = n_comps == 1;
    let norm_ok = n.inward_faces == 0;
    let degen_ok = n.degenerate_faces == 0;
    let any_issue = !wt.is_watertight || !chi_ok || !comps_ok || !norm_ok || !degen_ok
        || !positive || !below_naive;
    let genus = (2 - expected_chi) / 2;

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (naive sum {v_naive:.4})");
    println!(
        "    Bounds     : V > 0 [{}]  V < V_naive [{}]",
        if positive { "PASS" } else { "FAIL" },
        if below_naive { "PASS" } else { "WARN" },
    );
    println!(
        "    Watertight : {}  (boundary={}, non-manifold={})",
        wt.is_watertight, wt.boundary_edge_count, wt.non_manifold_edge_count
    );
    println!(
        "    Euler χ    : {:?}  (expected {expected_chi}, genus {genus})  [{}]",
        wt.euler_characteristic,
        if chi_ok { "PASS" } else { "WARN" }
    );
    println!(
        "    Components : {n_comps}  [{}]",
        if comps_ok { "PASS" } else { "WARN phantom islands" }
    );
    println!(
        "    Normals    : outward={}, inward={} ({:.1}%), degen={}  [{}]",
        n.outward_faces,
        n.inward_faces,
        if mesh.face_count() > 0 {
            n.inward_faces as Real / mesh.face_count() as Real * 100.0
        } else {
            0.0
        },
        n.degenerate_faces,
        if norm_ok && degen_ok { "PASS" } else { "WARN" }
    );
    println!(
        "    Alignment  : mean={:.4}  min={:.4}",
        n.face_vertex_alignment_mean, n.face_vertex_alignment_min
    );
    println!("    Elapsed    : {ms} ms");

    if any_issue {
        println!("    *** GEOMETRY ISSUES DETECTED ***");
        if !wt.is_watertight {
            println!(
                "       - Not watertight: {} boundary + {} non-manifold edge(s)",
                wt.boundary_edge_count, wt.non_manifold_edge_count
            );
        }
        if !chi_ok {
            println!(
                "       - Euler χ = {:?} (expected {expected_chi}): phantom islands or non-manifold topology",
                wt.euler_characteristic
            );
        }
        if !comps_ok {
            println!(
                "       - {} connected component(s): {} phantom island(s) present",
                n_comps, n_comps.saturating_sub(1)
            );
        }
        if !norm_ok {
            println!(
                "       - {}/{} face(s) with inward normals: winding order errors",
                n.inward_faces, mesh.face_count()
            );
        }
        if !degen_ok {
            println!(
                "       - {} degenerate face(s) with zero area",
                n.degenerate_faces
            );
        }
        if !positive {
            println!("       - Negative volume: mesh is inside-out or empty");
        }
    }
    println!();
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
