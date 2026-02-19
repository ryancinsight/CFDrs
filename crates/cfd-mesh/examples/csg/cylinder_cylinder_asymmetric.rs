//! CSG Cylinder–Cylinder (Asymmetric): union, intersection, and difference
//!
//! Two Y-axis cylinders with the **same radius** but **different heights** so
//! that their cap planes are at different Y coordinates — no coplanar boundary
//! faces exist.  This tests the standard Mesh Arrangement barrel-seam path
//! without any coplanar-cap complication.
//!
//! ## Geometry
//!
//! ```text
//! Cylinder A : base (−0.3, −1.5, 0), r = 0.6, h = 3.0 mm  →  Y ∈ [−1.5, 1.5]
//! Cylinder B : base (+0.3, −2.0, 0), r = 0.6, h = 4.0 mm  →  Y ∈ [−2.0, 2.0]
//!   axis separation d = 0.6 mm = r_A  →  θ = π/3
//!   caps at y = ±1.5 (A) and y = −2.0, +2.0 (B) — none coplanar
//! ```
//!
//! The overlap region is the cylindrical lens cross-section extruded over the
//! **shared height** Y ∈ [−1.5, 1.5] (A's full extent, since A is shorter):
//!
//! ```text
//! θ         = π / 3                (same as symmetric case, d = r)
//! A_segment = r²(θ − sin θ cos θ)
//! h_overlap = min(h_A_top, h_B_top) − max(h_A_bot, h_B_bot)
//!           = min(1.5, 2.0) − max(−1.5, −2.0) = 1.5 − (−1.5) = 3.0 mm
//! V_∩       = 2 · h_overlap · A_segment  ≈ 1.3267 mm³
//! ```
//!
//! | Operation | Expected (mm³)              | Pipeline    |
//! |-----------|----------------------------|-------------|
//! | A ∪ B     | V_A + V_B − V_∩            | Arrangement |
//! | A ∩ B     | V_∩ ≈ 1.3267               | Arrangement |
//! | A \ B     | V_A − V_∩ ≈ 2.0662         | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_asymmetric
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::io::stl;
use cfd_mesh::{IndexedMesh, analyze_normals};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (Asymmetric): Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — different heights, no coplanar caps)");
    println!("=================================================================");

    let r:    f64 = 0.6;   // both cylinder radii (equal)
    let h_a:  f64 = 3.0;   // A: Y ∈ [−1.5, 1.5]
    let h_b:  f64 = 4.0;   // B: Y ∈ [−2.0, 2.0]  — taller, caps offset from A
    let d:    f64 = r;     // axis separation = r  →  θ = π/3

    let v_a = std::f64::consts::PI * r * r * h_a;
    let v_b = std::f64::consts::PI * r * r * h_b;

    // Overlap height = intersection of Y extents:
    //   A: [−1.5, 1.5]   B: [−2.0, 2.0]   →  [−1.5, 1.5]
    let y_a_bot = -h_a / 2.0;  // −1.5
    let y_b_bot = -h_b / 2.0;  // −2.0
    let y_a_top =  h_a / 2.0;  // +1.5
    let y_b_top =  h_b / 2.0;  // +2.0
    let h_overlap = (y_a_top.min(y_b_top) - y_a_bot.max(y_b_bot)).max(0.0); // 3.0

    // Lens cross-section area (d = r → θ = π/3):
    let theta     = (d / (2.0 * r)).acos();                                  // π/3
    let a_seg     = r * r * (theta - theta.sin() * theta.cos());
    let v_intersect = 2.0 * h_overlap * a_seg;

    println!("  Cylinder A : base ({:.2},{:.2},0), r={r}, h={h_a}  V = {v_a:.4} mm³",
        -d / 2.0, y_a_bot);
    println!("  Cylinder B : base ({:.2},{:.2},0), r={r}, h={h_b}  V = {v_b:.4} mm³",
        d / 2.0, y_b_bot);
    println!("  Axis sep d={d} mm (=r → θ=π/3);  caps at y=±1.5 (A) vs y=±2.0 (B)");
    println!("  Overlap h={h_overlap} mm, lens cross-section V_∩ = {v_intersect:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", v_a + v_b - v_intersect);
    println!("    Intersection : {v_intersect:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_a - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    // A: base at (−d/2, y_a_bot, 0) so Y ∈ [−1.5, +1.5]
    let cyl_a = Cylinder {
        base_center: Point3r::new(-d / 2.0, y_a_bot, 0.0),
        radius: r,
        height: h_a,
        segments: 64,
    }.build()?;
    // B: base at (+d/2, y_b_bot, 0) so Y ∈ [−2.0, +2.0] — protrudes 0.5 mm past each A cap
    let cyl_b = Cylinder {
        base_center: Point3r::new(d / 2.0, y_b_bot, 0.0),
        radius: r,
        height: h_b,
        segments: 64,
    }.build()?;
    println!("  Mesh built: {} + {} cylinder faces  ({} ms)",
        cyl_a.face_count(), cyl_b.face_count(), t_build.elapsed().as_millis());
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Union (A ∪ B)", &mut result, v_a + v_b - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_asymmetric_union.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_asymmetric_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_asymmetric_intersection.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_asymmetric_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_a - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_asymmetric_difference.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_asymmetric_difference.stl");
        println!();
    }

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

