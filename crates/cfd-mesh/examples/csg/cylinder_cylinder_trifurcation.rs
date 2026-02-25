//! CSG Cylinder–Cylinder (Trifurcation): union and trunk difference
//!
//! One **trunk** cylinder along +X ends at a junction from which three
//! **branches** diverge in the XY plane:
//!
//! - Branch B (upper-forward): at +θ from the +X axis
//! - Branch C (perpendicular): at +90° from the +X axis
//! - Branch D (lower-forward): at −θ from the +X axis
//!
//! ```text
//!               ╱── B (+θ, upper-forward)
//!              │    C (+90°, straight up)
//!  ──────A────⊙
//!              ╲── D (−θ, lower-forward)
//!  ↑ inlet      ↑ junction (0, 0, 0)
//! ```
//!
//! At θ = 45° the three branches are evenly spaced 45° apart.
//! Sweeping θ from 30° → 60° narrows / widens the outer branches while the
//! perpendicular prong remains fixed.
//!
//! ## Degeneracy avoidance
//!
//! A "flat trident" with a centre-forward prong at θ = 0° places one branch
//! coaxial with the trunk.  Same-radius, same-axis tubes produce coincident
//! lateral surfaces and `csg_boolean_indexed(Union, …)` returns `NotWatertight`.
//! Choosing the perpendicular branch at +90° avoids this entirely.
//!
//! The CSG boolean kernel uses a per-face 2-D CDT to corefine intersection
//! curves.  All four cylinders lie in the XY plane, so face normals point
//! along ±Z and the CDT projections are always consistent — no 2-D segment
//! crossing arises.  (Three branches distributed in 3-D at 120° azimuthal
//! spacing would cause the CDT to fail with "segments intersect in their
//! interiors" because out-of-plane face projections can cross.)
//!
//! ## Operations shown
//!
//! | Operation                        | Variants           |
//! |----------------------------------|--------------------|
//! | Union `((A∪B)∪C)∪D`             | θ = 30°, 45°, 60° |
//! | Trunk difference `A \ (B∪C∪D)`  | θ = 45°            |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_trifurcation
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

/// Tube radius [mm]. Shared by trunk and all three branches.
const R: f64 = 0.5;

/// Trunk centreline length [mm] (from inlet to nominal junction).
const H_TRUNK: f64 = 3.0;

/// Branch length [mm] (from junction to branch outlet).
const H_BRANCH: f64 = 3.0;

/// Trunk extension past the junction [mm].
///
/// Extends the trunk by ε so that branches overlap the trunk barrel rather
/// than only touching at a single point. Small relative to the tube radius
/// (~10 % of R) so it introduces negligible volume error.
const EPS: f64 = R * 0.10; // 0.05 mm

/// Circumferential segments for all cylinders.
const SEGS: usize = 64;

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (Trifurcation)");
    println!("  Outer branches ±θ  |  Perpendicular centre branch +90°");
    println!("  Union | Trunk Difference");
    println!(
        "  r={R} mm  h_trunk={H_TRUNK} mm  h_branch={H_BRANCH} mm  ε={EPS:.3} mm"
    );
    println!("=================================================================");
    println!();

    let v_a = std::f64::consts::PI * R * R * (H_TRUNK + EPS);
    let v_b = std::f64::consts::PI * R * R * H_BRANCH;
    let v_naive = v_a + 3.0 * v_b;

    println!("  V_A (trunk+ε)       = πr²(h_trunk + ε) = {v_a:.4} mm³");
    println!("  V_B = V_C = V_D     = πr²·h_branch      = {v_b:.4} mm³");
    println!("  V_naive (sum)       = V_A + 3·V_B        = {v_naive:.4} mm³");
    println!("  Actual union        < V_naive (overlaps subtracted)");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // ── Union sweep at 30°, 45°, 60° ─────────────────────────────────────────

    for &theta_deg in &[30.0_f64, 45.0, 60.0] {
        let theta = theta_deg.to_radians();

        println!(
            "  ─── θ = {theta_deg}°  (outer branches ±{theta_deg}°, centre +90°) ─────────────"
        );

        let t_build = Instant::now();
        let (trunk, b_up, b_perp, b_dn) = build_trifurcation(theta)?;
        println!(
            "  Meshes built: {} + {} + {} + {} faces  ({} ms)",
            trunk.face_count(),
            b_up.face_count(),
            b_perp.face_count(),
            b_dn.face_count(),
            t_build.elapsed().as_millis()
        );

        // Union: ((A ∪ B) ∪ C) ∪ D
        let t0 = Instant::now();
        let a_union_b = csg_boolean_indexed(BooleanOp::Union, &trunk, &b_up)?;
        let ab_union_c = csg_boolean_indexed(BooleanOp::Union, &a_union_b, &b_perp)?;
        let mut result = csg_boolean_indexed(BooleanOp::Union, &ab_union_c, &b_dn)?;
        let ms = t0.elapsed().as_millis();

        report_union(
            &format!("Union (A ∪ B ∪ C ∪ D)  θ={theta_deg}°"),
            &mut result,
            v_naive,
            ms,
        );

        let stl_name =
            format!("cylinder_cylinder_trifurcation_union_{:.0}deg.stl", theta_deg);
        write_stl(&result, &out_dir.join(&stl_name))?;
        println!("  STL: outputs/csg/{stl_name}");
        println!();
    }

    // ── Trunk difference A \ (B ∪ C ∪ D) at θ = 45° ─────────────────────────

    println!("  ─── θ = 45°: trunk difference ───────────────────────────────────");
    println!();

    {
        let (trunk_45, b_up_45, b_perp_45, b_dn_45) =
            build_trifurcation(std::f64::consts::FRAC_PI_4)?;

        let t0 = Instant::now();
        let b_up_perp =
            csg_boolean_indexed(BooleanOp::Union, &b_up_45, &b_perp_45)?;
        let branches = csg_boolean_indexed(BooleanOp::Union, &b_up_perp, &b_dn_45)?;
        let mut result =
            csg_boolean_indexed(BooleanOp::Difference, &trunk_45, &branches)?;
        let ms = t0.elapsed().as_millis();

        // Three branches clip slightly more junction material than the 2-branch
        // y-junction, so tolerance is 20 % (vs 15 % for the y-junction trunk diff).
        report(
            "Trunk Difference (A \\ (B ∪ C ∪ D))  θ=45°",
            &mut result,
            v_a,
            0.20,
            ms,
        );

        write_stl(
            &result,
            &out_dir
                .join("cylinder_cylinder_trifurcation_trunk_diff_45deg.stl"),
        )?;
        println!(
            "  STL: outputs/csg/cylinder_cylinder_trifurcation_trunk_diff_45deg.stl"
        );
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Trifurcation builder ───────────────────────────────────────────────────────

/// Build the four component cylinders for a planar trifurcation with outer
/// branch half-angle `theta` (radians, measured from the +X trunk axis).
///
/// Returns `(trunk A, branch_up B, branch_perp C, branch_dn D)`.
///
/// # Layout (XY plane)
///
/// ```text
///               ╱── B at +θ
///              │    C at +90° (perpendicular)
///  ────A──────⊙
///              ╲── D at −θ
/// ```
///
/// All cylinders are built in canonical `+Y` orientation, then redirected by
/// rotating about Z:
///
/// | Part       | Rotation about Z    | Translation           | World axis         |
/// |------------|---------------------|-----------------------|--------------------|
/// | Trunk A    | −π/2 (+Y → +X)     | (−H_TRUNK, 0, 0)      | +X                 |
/// | Branch B   | θ − π/2            | (0, 0, 0)             | (cosθ, sinθ, 0)    |
/// | Branch C   | 0  (+Y stays +Y)   | (0, 0, 0)             | +Y                 |
/// | Branch D   | −θ − π/2           | (0, 0, 0)             | (cosθ, −sinθ, 0)   |
///
/// The trunk is extended by `EPS` past the junction so each branch cylinder
/// has genuine material overlap with the trunk barrel.
fn build_trifurcation(
    theta: f64,
) -> Result<
    (IndexedMesh, IndexedMesh, IndexedMesh, IndexedMesh),
    Box<dyn std::error::Error>,
> {
    // ── Trunk A: (−H_TRUNK, 0, 0) → (EPS, 0, 0) along +X ────────────────────
    let trunk = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: H_TRUNK + EPS,
            segments: SEGS,
        }
        .build()?;
        // Rotate −π/2 about Z: +Y → +X.
        // After rotation the cylinder spans (0,0,0) → (H_TRUNK+EPS, 0, 0).
        // Translate (−H_TRUNK, 0, 0) so the inlet is at (−H_TRUNK, 0, 0)
        // and the extended tip is at (EPS, 0, 0).
        let rot =
            UnitQuaternion::<Real>::from_axis_angle(&Vector3::z_axis(), -FRAC_PI_2);
        let iso =
            Isometry3::from_parts(Translation3::new(-H_TRUNK, 0.0, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }
        .evaluate()?
    };

    // ── Branches: all built as canonical +Y cylinders, then rotated ───────────
    //
    // A +Y cylinder rotated by angle α about Z ends up pointing along the
    // direction (cos(α + π/2), sin(α + π/2), 0) = (−sinα, cosα, 0).
    // Equivalently: rotating by (φ − π/2) points the cylinder at angle φ from +X.

    // Branch B: at +θ from +X
    let b_up = make_branch_planar(theta)?;

    // Branch C: at +90° from +X (straight up along +Y).
    // Rotation = 90° − π/2 = 0 → identity; the canonical +Y cylinder is already correct.
    let b_perp = make_branch_planar(FRAC_PI_2)?;

    // Branch D: at −θ from +X (mirror of B below the trunk axis)
    let b_dn = make_branch_planar(-theta)?;

    Ok((trunk, b_up, b_perp, b_dn))
}

/// Build a single branch cylinder at `angle_from_x` (radians) from the +X axis
/// in the XY plane, base at the junction origin.
///
/// The canonical `+Y` cylinder is redirected by rotating `(angle_from_x − π/2)`
/// about the Z axis, which maps +Y to `(cos angle_from_x, sin angle_from_x, 0)`.
fn make_branch_planar(
    angle_from_x: f64,
) -> Result<IndexedMesh, Box<dyn std::error::Error>> {
    let raw = Cylinder {
        base_center: Point3r::new(0.0, 0.0, 0.0),
        radius: R,
        height: H_BRANCH,
        segments: SEGS,
    }
    .build()?;
    let rot = UnitQuaternion::<Real>::from_axis_angle(
        &Vector3::z_axis(),
        angle_from_x - FRAC_PI_2,
    );
    let iso = Isometry3::from_parts(Translation3::new(0.0, 0.0, 0.0), rot);
    Ok(CsgNode::Transform {
        node: Box::new(CsgNode::Leaf(raw)),
        iso,
    }
    .evaluate()?)
}

// ── Report helpers ─────────────────────────────────────────────────────────────

/// Report for union operations: checks bounds `0 < V < V_naive` instead of a
/// specific expected value (no simple closed form for 3-branch overlap).
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

/// Report for operations with an approximate expected volume (difference).
/// Volume error is compared against `expected` with tolerance `tol` (fraction,
/// e.g. 0.20 = 20 %).
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
/// For N genus-0 closed components: χ = V − E + F = 2N, so χ > 2 reveals
/// phantom islands that `is_watertight` alone cannot catch.
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
