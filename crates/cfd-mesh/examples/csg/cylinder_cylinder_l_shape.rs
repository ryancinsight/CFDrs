//! CSG / Primitive Cylinder–Cylinder (L-shape): sharp corner vs rounded corner
//!
//! Two equal-radius cylinders arranged as an **L**: one runs along +Y (the
//! "stem"), the other runs along +X starting from the top of the stem (the
//! "arm").  Their axes meet at a right angle at the shared endpoint.
//!
//! The example demonstrates **both** construction strategies, controlled by the
//! `ROUNDED` constant:
//!
//! ## Sharp corner (`ROUNDED = false`)
//!
//! The two cylinders are unioned with CSG.  Their end-caps intersect at the
//! corner — the top cap of the stem and the base cap of the arm share the same
//! quarter-circle footprint.  The resulting solid has a sharp internal corner
//! where the two cylinder walls meet.
//!
//! ```text
//!   ┌──────────────────────────┐   y = h
//!   │   arm (axis +X)          │
//!   └──────────────────────────┘
//!   │
//!   │  stem (axis +Y)
//!   │
//!   └   y = 0
//! ```
//!
//! Geometry:
//! ```text
//! Stem (A): base (0, 0, 0), axis +Y, r, h
//! Arm  (B): base (0, h, 0), axis +X, r, h   (B is A rotated −90°Z then translated (0,h,0))
//! V_∩ = (1/4) · (16r³/3) = 4r³/3           (quarter-Steinmetz at the inside corner)
//! ```
//!
//! ## Rounded corner (`ROUNDED = true`)
//!
//! Instead of CSG, the L-shape is built from **three watertight primitives**
//! unioned together:
//!
//! 1. **Stem straight** — a cylinder running from `(0, 0, 0)` to `(0, h−R, 0)`,
//!    stopping short of the corner by one bend-radius `R`.
//! 2. **Corner elbow** — a 90° `Elbow` with `bend_radius = R`, sweeping from the
//!    +Y direction to the +X direction.  Its inlet sits at `(0, h−R, 0)` and its
//!    outlet at `(R, h, 0)` (after coordinate transform).  The elbow lives in the
//!    XY plane; the Elbow primitive lives in the XZ plane, so a 90° rotation about
//!    Y is applied.
//! 3. **Arm straight** — a cylinder from `(R, h, 0)` to `(h, h, 0)` (the
//!    remaining arm length after the bend).
//!
//! The three pieces are unioned with CSG.  Because the cylinder end-caps sit
//! flush against the elbow's inlet/outlet caps, the union is a clean merge with
//! no material overlap — only the shared cap faces need to be resolved.
//!
//! ```text
//!           ╭──────────────────╮   y = h
//!          ╱  arm cylinder      ╲
//!         │                      │
//!          ╲                    ╱
//!       ╭──╯  elbow corner  ╰──╮  y = h-R
//!       │                       │
//!       │  stem cylinder        │
//!       │                       │
//!       ╰───────────────────────╯  y = 0
//! ```
//!
//! ## Why the rounded approach is better for CFD
//!
//! Sharp internal corners create singular pressure and velocity fields in
//! Navier–Stokes solvers.  A rounded corner with bend-radius R ≥ r avoids
//! the singularity and gives a mesh-quality-friendly transition.  The elbow
//! is the canonical building block for millifluidic channel networks.
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_l_shape
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

use cfd_mesh::core::scalar::{Point3r, Real};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::csg::CsgNode;
use cfd_mesh::geometry::primitives::{Cylinder, Elbow, PrimitiveMesh};
use cfd_mesh::io::stl;
use cfd_mesh::{IndexedMesh, analyze_normals};

// ── Configuration ─────────────────────────────────────────────────────────────

/// Switch between the two construction strategies.
///
/// `false` — CSG union of two cylinders with a sharp internal corner.
/// `true`  — Three-piece assembly (stem + elbow + arm) with a smooth rounded corner.
const ROUNDED: bool = true;

// ── Geometry parameters ────────────────────────────────────────────────────────

/// Tube radius [mm].  Shared by stem, arm, and elbow.
const R: f64 = 0.5;

/// Length of each straight arm [mm] (stem height, arm length).
const H: f64 = 3.0;

/// Bend-centreline radius for the rounded corner [mm].
/// Must be > R (tube radius).  R_BEND = 2·R gives a moderate fillet.
const R_BEND: f64 = 2.0 * R;  // = 1.0 mm

// ── Main ──────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    if ROUNDED {
        run_rounded()
    } else {
        run_sharp()
    }
}

// ── Sharp corner (CSG union) ──────────────────────────────────────────────────

fn run_sharp() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  L-shape: SHARP corner (CSG union of two cylinders)");
    println!("=================================================================");

    let v_cyl = std::f64::consts::PI * R * R * H;
    // The intersection at the inside corner is 1/4 of the Steinmetz bicylinder.
    // Both cylinders have radius R; their axes meet at a right angle at one end.
    // The overlapping volume is the quarter of both cylinders that lies within
    // the footprint of the other: V_∩ = (1/4)(16R³/3) = 4R³/3.
    let v_intersect = 4.0 / 3.0 * R * R * R;

    println!("  Stem (A): base (0,0,0), axis +Y, r={R}, h={H}  V = {v_cyl:.4} mm³");
    println!("  Arm  (B): base (0,{H},0), axis +X, r={R}, h={H}  V = {v_cyl:.4} mm³");
    println!("  Corner overlap (quarter-Steinmetz): V_∩ = 4r³/3 = {v_intersect:.4} mm³");
    println!();
    println!("  Expected union volume: {:.4} mm³", 2.0 * v_cyl - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t0 = Instant::now();

    // Stem: Y-axis, base at origin.
    let stem = Cylinder {
        base_center: Point3r::new(0.0, 0.0, 0.0),
        radius: R,
        height: H,
        segments: 64,
    }.build()?;

    // Arm: Y-axis cylinder rotated −90° about Z (+Y→+X), then translated so
    // its base (which becomes its −X end after rotation) sits at (0, H, 0).
    // After −90°Z rotation the cylinder spans X ∈ [0, H] at y = H, z = 0.
    let arm = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: H,
            segments: 64,
        }.build()?;
        // Rotate −90° about Z: +Y → +X.
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,
        );
        // Translate to corner: (0, H, 0).
        let iso = Isometry3::from_parts(Translation3::new(0.0, H, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }.evaluate()?
    };

    let build_ms = t0.elapsed().as_millis();
    println!("  Meshes built: {} + {} faces  ({build_ms} ms)", stem.face_count(), arm.face_count());
    println!();

    // Union
    let t1 = Instant::now();
    let mut result = csg_boolean_indexed(BooleanOp::Union, &stem, &arm)?;
    let union_ms = t1.elapsed().as_millis();
    report("Union (stem ∪ arm) — sharp corner", &mut result, 2.0 * v_cyl - v_intersect, 0.05, union_ms);
    let stl_path = out_dir.join("cylinder_cylinder_l_shape_sharp.stl");
    write_stl(&result, &stl_path)?;
    println!("  STL: outputs/csg/cylinder_cylinder_l_shape_sharp.stl");

    println!("=================================================================");
    Ok(())
}

// ── Rounded corner (stem + elbow + arm) ──────────────────────────────────────

fn run_rounded() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  L-shape: ROUNDED corner (stem cylinder + 90° elbow + arm cylinder)");
    println!("=================================================================");

    // The straight legs are shortened by R_BEND so the elbow fits at the corner.
    let straight_len = H - R_BEND;

    // Each piece slightly overlaps its neighbour so that the CSG union has
    // genuine intersecting geometry at each joint rather than exactly coplanar
    // coincident caps (which the arrangement pipeline cannot cleanly resolve).
    // EPS is small relative to the tube radius so the volume error it introduces
    // is negligible (~EPS / straight_len < 0.1%).
    let eps = R * 0.05;  // 5% of tube radius ≈ 0.025 mm

    // Expected volumes (additive, overlaps ≈ π r² eps × 2 junctions → negligible):
    //   Stem:  π r² · (straight_len + eps)
    //   Elbow: π r² · R_BEND · (π/2)   — Pappus, arc length only
    //   Arm:   π r² · (straight_len + eps)
    //   Total ≈ 2·π r²·straight_len + π r²·R_BEND·(π/2)  (eps terms cancel in union)
    let v_straight = std::f64::consts::PI * R * R * straight_len;
    let v_elbow    = std::f64::consts::PI * R * R * R_BEND * std::f64::consts::FRAC_PI_2;
    let v_total    = 2.0 * v_straight + v_elbow;

    println!("  Tube radius r = {R} mm,  bend radius R_bend = {R_BEND} mm");
    println!("  Arm length H = {H} mm  →  straight leg = {straight_len} mm (+{eps:.3} mm overlap)");
    println!();
    println!("  Piece volumes (Pappus, ignoring eps overlaps):");
    println!("    Stem straight  : {v_straight:.4} mm³");
    println!("    Arm  straight  : {v_straight:.4} mm³");
    println!("    90° elbow      : {v_elbow:.4} mm³  (π r² R_bend π/2)");
    println!("    Total (approx) : {v_total:.4} mm³");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t0 = Instant::now();

    // ── Piece 1: Stem straight ────────────────────────────────────────────────
    // Y-axis, from (0, 0, 0) up to (0, straight_len + eps, 0).
    // The extra eps pushes the stem's top cap just inside the elbow barrel,
    // giving the CSG union a real intersecting surface to cut.
    let stem = Cylinder {
        base_center: Point3r::new(0.0, 0.0, 0.0),
        radius: R,
        height: straight_len + eps,
        segments: 64,
    }.build()?;

    // ── Piece 2: Elbow (90° bend) ─────────────────────────────────────────────
    //
    // Elbow primitive in XZ:   C(α) = (R_BEND·(1−cosα), 0, R_BEND·sinα)
    //   inlet  at (0, 0, 0),            tangent T(0)   = +Z
    //   outlet at (R_BEND, 0, R_BEND),  tangent T(π/2) = +X
    //
    // Target orientation: inlet tangent = +Y, outlet tangent = +X.
    // Rotation −90° about X maps (x,y,z) → (x, z, −y):
    //   +Z → +Y  ✓     +X → +X  ✓
    //   outlet (R_BEND,0,R_BEND) → (R_BEND, R_BEND, 0)
    //
    // Translate inlet to (0, straight_len, 0) — the elbow begins at the
    // nominal (non-overlapped) stem top so its body overlaps the stem top:
    //   inlet  → (0,        straight_len,         0)
    //   outlet → (R_BEND,   R_BEND + straight_len, 0) = (R_BEND, H, 0)
    //
    // Final world geometry:
    //   Stem:  (0,0,0) → (0, straight_len+eps, 0)        along +Y
    //   Elbow: arc from (0, straight_len, 0) → (R_BEND, H, 0)   in XY
    //   Arm:   (R_BEND−eps, H, 0) → (R_BEND + straight_len, H, 0)  along +X
    let arm_y = H;   // y-level of the arm centreline

    let elbow_mesh = {
        let raw = Elbow {
            tube_radius:   R,
            bend_radius:   R_BEND,
            bend_angle:    std::f64::consts::FRAC_PI_2,
            tube_segments: 64,
            arc_segments:  32,
        }.build()?;
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::x_axis(),
            -std::f64::consts::FRAC_PI_2,  // +Z → +Y, +X → +X
        );
        let iso = Isometry3::from_parts(Translation3::new(0.0, straight_len, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }.evaluate()?
    };

    // ── Piece 3: Arm straight ─────────────────────────────────────────────────
    // X-axis, from (R_BEND − eps, arm_y, 0) along +X for (straight_len + eps) mm.
    // Starting eps before the elbow outlet so the arm's base cap intersects the
    // elbow barrel rather than sitting flush against it.
    let arm = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, 0.0, 0.0),
            radius: R,
            height: straight_len + eps,
            segments: 64,
        }.build()?;
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,  // +Y → +X
        );
        // Start eps before the elbow outlet so arm base pierces the elbow barrel.
        let iso = Isometry3::from_parts(Translation3::new(R_BEND - eps, arm_y, 0.0), rot);
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }.evaluate()?
    };

    let build_ms = t0.elapsed().as_millis();
    println!("  Pieces built: stem={} elbow={} arm={} faces  ({build_ms} ms)",
        stem.face_count(), elbow_mesh.face_count(), arm.face_count());
    println!();

    // ── Union the three pieces ────────────────────────────────────────────────
    // Each pair overlaps by eps, so the union correctly trims the interior caps.
    // Sequence: (stem ∪ elbow) ∪ arm.
    let t1 = Instant::now();
    let stem_elbow = csg_boolean_indexed(BooleanOp::Union, &stem, &elbow_mesh)?;
    let mut result = csg_boolean_indexed(BooleanOp::Union, &stem_elbow, &arm)?;
    let union_ms = t1.elapsed().as_millis();

    // Tolerance on expected: the eps overlaps add ~2·π r² eps of material that
    // cancels in the union, so the true volume is still ≈ v_total. Use 5% tol.
    report("Union (stem ∪ elbow ∪ arm) — rounded corner", &mut result, v_total, 0.05, union_ms);
    let stl_path = out_dir.join("cylinder_cylinder_l_shape_rounded.stl");
    write_stl(&result, &stl_path)?;
    println!("  STL: outputs/csg/cylinder_cylinder_l_shape_rounded.stl");

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol    = mesh.signed_volume();
    let is_wt  = mesh.is_watertight();
    let n      = analyze_normals(mesh);
    let err    = (vol - expected).abs() / expected.abs().max(1e-12);
    let status = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected {expected:.4})");
    println!("    Vol error  : {:.2}%  [{status}]", err * 100.0);
    println!("    Watertight : {is_wt}");
    println!("    Normals    : outward={}, inward={} ({:.1}%), degen={}",
        n.outward_faces, n.inward_faces,
        if mesh.face_count() > 0 { n.inward_faces as Real / mesh.face_count() as Real * 100.0 } else { 0.0 },
        n.degenerate_faces);
    println!("    Alignment  : mean={:.4}  min={:.4}",
        n.face_vertex_alignment_mean, n.face_vertex_alignment_min);
    println!("    Elapsed    : {ms} ms");
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file  = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
