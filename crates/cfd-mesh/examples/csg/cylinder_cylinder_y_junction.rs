//! CSG Cylinder–Cylinder (Y-Junction): union, branch intersection, trunk difference
//!
//! One **trunk** cylinder along +X ends at a **junction point** from which two
//! symmetric **branches** diverge at half-angle θ above and below the +X axis.
//! The trunk is extended by a small ε past the junction so that each branch
//! cylinder has genuine material overlap with the trunk, avoiding degenerate
//! tangent-only touching.
//!
//! ## Geometry (example at θ = 45°)
//!
//! ```text
//!                    ╱── branch_up ──▶ (cos θ · L, +sin θ · L, 0)
//!  ──────trunk──── ⊙
//!  (-L, 0, 0)       ╲── branch_dn ──▶ (cos θ · L, −sin θ · L, 0)
//!  ↑ inlet           ↑ junction (0, 0, 0)
//! ```
//!
//! - Trunk (A): axis +X, from `(−h_trunk, 0, 0)` to `(ε, 0, 0)`.
//!   Built as a +Y cylinder, rotated −90° about Z (+Y → +X), translated (−h_trunk, 0, 0).
//! - Branch up (B): axis `(cos θ, +sin θ, 0)`, base at origin, length h_branch.
//!   Built as +Y cylinder, rotated `(θ − π/2)` about Z.
//! - Branch down (C): axis `(cos θ, −sin θ, 0)`, base at origin, length h_branch.
//!   Built as +Y cylinder, rotated `(−θ − π/2)` about Z.
//!
//! ## Volume analysis
//!
//! Let V_naive = πr²(h_trunk + ε) + 2πr²·h_branch (sum ignoring intersections).
//!
//! The union subtracts the branch–branch overlap (V_B∩C) and the small
//! trunk–branch overlaps (≈ ε-scale), so V_union < V_naive.
//!
//! Both branches are semi-infinite (start at origin, extend outward), so B ∩ C
//! is the quarter of the Steinmetz bicylinder where both axial projections ≥ 0:
//!   θ = 45° (90° between axes): V_B∩C = 4r³/3  (quarter-Steinmetz)
//!
//! ## Operations shown
//!
//! | Operation          | Meshes           | Variants        |
//! |--------------------|------------------|-----------------|
//! | Union A ∪ B ∪ C   | trunk + both branches | θ = 30°, 45°, 60° |
//! | Intersection B ∩ C | branches only    | θ = 45°         |
//! | Difference A \ BC  | trunk minus both | θ = 45°         |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_y_junction
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::f64::consts::FRAC_PI_2;
use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::application::csg::CsgNode;
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::domain::topology::adjacency::AdjacencyGraph;
use cfd_mesh::domain::topology::connectivity::connected_components;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, IndexedMesh};

// ── Geometry parameters ────────────────────────────────────────────────────────

/// Tube radius [mm].  Shared by trunk and both branches.
const R: f64 = 0.5;

/// Trunk centreline length [mm] (from inlet to nominal junction).
const H_TRUNK: f64 = 3.0;

/// Branch length [mm] (from junction to branch outlet).
const H_BRANCH: f64 = 3.0;

/// Trunk extension past the junction [mm].
///
/// Extends the trunk by ε so that branches overlap the trunk barrel rather
/// than only touching at a single point.  Small relative to the tube radius
/// (~10 % of R) so it introduces negligible volume error.
const EPS: f64 = R * 0.10; // 0.05 mm

/// Circumferential segments for all cylinders.
const SEGS: usize = 64;

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (Y-Junction)");
    println!("  Union | Branch Intersection | Trunk Difference");
    println!("  r={R} mm  h_trunk={H_TRUNK} mm  h_branch={H_BRANCH} mm  ε={EPS:.3} mm");
    println!("=================================================================");
    println!();

    let v_a = std::f64::consts::PI * R * R * (H_TRUNK + EPS);
    let v_b = std::f64::consts::PI * R * R * H_BRANCH;
    let v_naive = v_a + 2.0 * v_b;

    println!("  V_A (trunk+ε)  = πr²(h_trunk + ε) = {v_a:.4} mm³");
    println!("  V_B = V_C      = πr²·h_branch      = {v_b:.4} mm³");
    println!("  V_naive (sum)  = V_A + 2·V_B       = {v_naive:.4} mm³");
    println!("  Actual union   < V_naive (overlap subtracted)");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // ── Union sweep at 30°, 45°, 60° ─────────────────────────────────────────

    for &theta_deg in &[30.0_f64, 45.0, 60.0] {
        let theta = theta_deg.to_radians();

        println!(
            "  ─── θ = {theta_deg}°  (branches {:.0}° apart) ─────────────────────────",
            theta_deg * 2.0
        );

        let t_build = Instant::now();
        let (trunk, branch_up, branch_dn) = build_y_junction(theta)?;
        println!(
            "  Meshes built: {} + {} + {} faces  ({} ms)",
            trunk.face_count(),
            branch_up.face_count(),
            branch_dn.face_count(),
            t_build.elapsed().as_millis()
        );

        // Union: (A ∪ B) ∪ C
        let t0 = Instant::now();
        let a_union_b = csg_boolean_indexed(BooleanOp::Union, &trunk, &branch_up)?;
        let mut result = csg_boolean_indexed(BooleanOp::Union, &a_union_b, &branch_dn)?;
        let ms = t0.elapsed().as_millis();

        report_union(
            &format!("Union (A ∪ B ∪ C)  θ={theta_deg}°"),
            &mut result,
            v_naive,
            ms,
        );

        let stl_name = format!("cylinder_cylinder_y_junction_union_{theta_deg:.0}deg.stl");
        write_stl(&result, &out_dir.join(&stl_name))?;
        println!("  STL: outputs/csg/{stl_name}");
        println!();
    }

    // ── 45° detailed operations ───────────────────────────────────────────────

    println!("  ─── θ = 45°: branch intersection and trunk difference ────────────");
    println!();

    let (trunk_45, branch_up_45, branch_dn_45) =
        build_y_junction(std::f64::consts::FRAC_PI_4)?;

    // ── Intersection B ∩ C ────────────────────────────────────────────────────
    // Both branches are semi-infinite cylinders starting at origin, axes at 90°.
    // The intersection is the quarter of the Steinmetz bicylinder where both
    // axial projections are non-negative:
    //
    //   V_full_Steinmetz  = 16r³/3  (two infinite perpendicular equal cylinders)
    //   V_B∩C             = (1/4)·(16r³/3) = 4r³/3  (quarter-Steinmetz)
    //
    // This matches the L-shape quarter-Steinmetz formula at a right-angle corner.
    {
        let v_quarter_steinmetz = 4.0 * R * R * R / 3.0;
        println!("  Quarter-Steinmetz (B ∩ C at 90°): V = 4r³/3 = {v_quarter_steinmetz:.4} mm³");

        let t0 = Instant::now();
        let mut result =
            csg_boolean_indexed(BooleanOp::Intersection, &branch_up_45, &branch_dn_45)?;
        let ms = t0.elapsed().as_millis();

        report(
            "Branch Intersection (B ∩ C)  θ=45°",
            &mut result,
            v_quarter_steinmetz,
            0.05, // tight tolerance: quarter-Steinmetz formula is accurate here
            ms,
        );

        write_stl(
            &result,
            &out_dir.join("cylinder_cylinder_y_junction_branch_intersection_45deg.stl"),
        )?;
        println!(
            "  STL: outputs/csg/cylinder_cylinder_y_junction_branch_intersection_45deg.stl"
        );
        println!();
    }

    // ── Difference A \ (B ∪ C) ────────────────────────────────────────────────
    // Removes both branches from the trunk.  Branches only overlap the trunk
    // in the ε-extension region (x ∈ [0, ε]), so V_diff ≈ V_A with small
    // reduction from the junction clips.
    {
        let t0_union = Instant::now();
        let branches_union =
            csg_boolean_indexed(BooleanOp::Union, &branch_up_45, &branch_dn_45)?;
        let t0_diff = Instant::now();
        let mut result =
            csg_boolean_indexed(BooleanOp::Difference, &trunk_45, &branches_union)?;
        let ms = t0_diff.elapsed().as_millis() + t0_union.elapsed().as_millis();

        report(
            "Trunk Difference (A \\ (B ∪ C))  θ=45°",
            &mut result,
            v_a,
            0.15, // branches clip only ε-scale material from trunk tip
            ms,
        );

        write_stl(
            &result,
            &out_dir.join("cylinder_cylinder_y_junction_trunk_diff_45deg.stl"),
        )?;
        println!("  STL: outputs/csg/cylinder_cylinder_y_junction_trunk_diff_45deg.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Y-junction builder ────────────────────────────────────────────────────────

/// Build the three component cylinders for a Y-junction with branch half-angle
/// `theta` (radians, measured from the +X trunk axis).
///
/// Returns `(trunk A, branch_up B, branch_down C)`.
///
/// # Transforms
///
/// All cylinders are built in canonical +Y orientation then rotated into place:
///
/// | Part         | Rotation about Z   | Translation      | World span                    |
/// |-------------|--------------------|--------------------|-------------------------------|
/// | Trunk A     | −π/2 (+Y → +X)    | (−h_trunk, 0, 0)  | (−h_trunk, 0, 0) → (ε, 0, 0) |
/// | Branch up B | θ − π/2           | (0, 0, 0)          | (0,0,0) → (L·cosθ, L·sinθ, 0)|
/// | Branch dn C | −θ − π/2          | (0, 0, 0)          | (0,0,0) → (L·cosθ,−L·sinθ, 0)|
fn build_y_junction(
    theta: f64,
) -> Result<(IndexedMesh, IndexedMesh, IndexedMesh), Box<dyn std::error::Error>> {
    // ── Trunk A: from (−H_TRUNK, 0, 0) to (EPS, 0, 0) along +X ──────────────
    let trunk = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: H_TRUNK + EPS,
            segments: SEGS,
        }
        .build()?;
        // Rotate −π/2 about Z: +Y → +X.  After rotation, cylinder spans
        // (0,0,0) → (H_TRUNK+EPS, 0, 0) along +X.
        // Translate (−H_TRUNK, 0, 0) so the inlet end is at (−H_TRUNK, 0, 0)
        // and the extended tip is at (EPS, 0, 0).
        let rot =
            UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -FRAC_PI_2);
        let iso = Isometry3::from_parts(Translation3::new(-H_TRUNK, 0.0, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }
        .evaluate()?
    };

    // ── Branch up B: from origin at angle +θ from +X ─────────────────────────
    let branch_up = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: H_BRANCH,
            segments: SEGS,
        }
        .build()?;
        // Rotate (θ − π/2) about Z.  At θ = π/2, no rotation (+Y stays +Y).
        // At θ = 0, rotation = −π/2 (+Y → +X, collinear with trunk — avoid!).
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            theta - FRAC_PI_2,
        );
        let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }
        .evaluate()?
    };

    // ── Branch down C: from origin at angle −θ from +X ───────────────────────
    let branch_dn = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: H_BRANCH,
            segments: SEGS,
        }
        .build()?;
        // Rotate (−θ − π/2) about Z.
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -theta - FRAC_PI_2,
        );
        let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }
        .evaluate()?
    };

    Ok((trunk, branch_up, branch_dn))
}

// ── Report helpers ────────────────────────────────────────────────────────────

/// Report for union operations: checks bounds (0 < V < V_naive) instead of a
/// specific expected value, since the branch–branch overlap has no simple
/// closed form.
fn report_union(label: &str, mesh: &mut IndexedMesh, v_naive: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let is_wt = mesh.is_watertight();
    let n = analyze_normals(mesh);

    let below_naive = vol < v_naive;
    let positive = vol > 0.0;
    let overlap_pct = (1.0 - vol / v_naive) * 100.0;

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!(
        "    Volume     : {vol:.4} mm³  (naive {v_naive:.4} mm³, overlap {overlap_pct:.1}%)"
    );
    println!(
        "    Bounds     : V > 0 [{}]  V < V_naive [{}]",
        pass(positive),
        pass(below_naive)
    );
    println!("    Watertight : {is_wt}");
    println!(
        "    Normals    : outward={}, inward={} ({:.1}%), degen={}",
        n.outward_faces,
        n.inward_faces,
        if mesh.face_count() > 0 {
            n.inward_faces as Real / mesh.face_count() as Real * 100.0
        } else {
            0.0
        },
        n.degenerate_faces
    );
    println!(
        "    Alignment  : mean={:.4}  min={:.4}",
        n.face_vertex_alignment_mean, n.face_vertex_alignment_min
    );
    println!("    Elapsed    : {ms} ms");
    connectivity_report(label, mesh, 1);
    println!();

    assert!(is_wt, "{label}: mesh must be watertight");
    assert!(positive, "{label}: mesh must have positive volume");
    assert!(below_naive, "{label}: union volume must be less than naive sum");
    assert_eq!(
        n.degenerate_faces, 0,
        "{label}: mesh must have no degenerate faces"
    );
    assert_eq!(
        n.inward_faces, 0,
        "{label}: mesh must have no inward-facing normals"
    );
}

/// Report for operations with an approximate expected volume (intersection,
/// difference).  Volume error is compared against `expected` with tolerance
/// `tol` (fraction, e.g. 0.15 = 15 %).
fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let is_wt = mesh.is_watertight();
    let n = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let status = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected ≈ {expected:.4})");
    println!("    Vol error  : {:.2}%  [{status}]", err * 100.0);
    println!("    Watertight : {is_wt}");
    println!(
        "    Normals    : outward={}, inward={} ({:.1}%), degen={}",
        n.outward_faces,
        n.inward_faces,
        if mesh.face_count() > 0 {
            n.inward_faces as Real / mesh.face_count() as Real * 100.0
        } else {
            0.0
        },
        n.degenerate_faces
    );
    println!(
        "    Alignment  : mean={:.4}  min={:.4}",
        n.face_vertex_alignment_mean, n.face_vertex_alignment_min
    );
    println!("    Elapsed    : {ms} ms");
    connectivity_report(label, mesh, 1);
    println!();

    assert!(is_wt, "{label}: mesh must be watertight");
    assert!(vol > 0.0, "{label}: mesh must have positive volume");
    assert_eq!(
        n.degenerate_faces, 0,
        "{label}: mesh must have no degenerate faces"
    );
    assert_eq!(
        n.inward_faces, 0,
        "{label}: orient_outward must eliminate all inward-facing normals"
    );
}

fn pass(ok: bool) -> &'static str {
    if ok { "PASS" } else { "FAIL" }
}

/// Print connected-component and Euler-characteristic diagnostics.
///
/// Detects isolated face islands that pass the watertight manifold check
/// because they are locally closed but disconnected from the main body.
/// For N genus-0 closed components: χ = V - E + F = 2N, so χ > 2 reveals
/// phantom islands that the `is_watertight` bool alone cannot catch.
fn connectivity_report(label: &str, mesh: &mut IndexedMesh, expected_components: usize) {
    mesh.rebuild_edges();
    let edges = mesh.edges_ref().unwrap();
    let wt = check_watertight(&mesh.vertices, &mesh.faces, edges);
    let adj = AdjacencyGraph::build(&mesh.faces, edges);
    let components = connected_components(&mesh.faces, &adj);

    let euler = wt.euler_characteristic.unwrap_or(i64::MIN);
    let expected_euler = 2 * expected_components as i64;

    println!("  ── Connectivity [{label}] ──");
    println!(
        "    Euler χ    : {} (expected {} for {} genus-0 body/bodies)",
        euler, expected_euler, expected_components,
    );
    println!("    Components : {}", components.len());
    for (i, comp) in components.iter().enumerate() {
        println!("      [{i}] {} faces", comp.len());
    }

    assert_eq!(
        components.len(),
        expected_components,
        "{label}: expected {expected_components} connected component(s), \
         got {} — retain_largest_component did not remove phantom islands",
        components.len(),
    );
    assert_eq!(
        wt.euler_characteristic,
        Some(expected_euler),
        "{label}: expected Euler χ = {expected_euler}, got {:?}",
        wt.euler_characteristic,
    );
}

fn write_stl(
    mesh: &IndexedMesh,
    path: &std::path::Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
