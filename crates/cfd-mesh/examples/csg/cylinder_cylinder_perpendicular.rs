//! CSG Cylinder–Cylinder (Perpendicular / Cross): union, intersection, and difference
//!
//! Two **equal-radius** cylinders whose axes are **perpendicular** and cross at
//! their mutual centres.  Cylinder A runs along **+Y**; Cylinder B runs along
//! **+X** (built as a Y-axis cylinder then rotated 90° about Z).
//!
//! ## Geometry
//!
//! ```text
//! Cylinder A : base (0, −1.5, 0),  r = 0.5 mm,  h = 3 mm  →  axis +Y
//! Cylinder B : base (0, −1.5, 0),  r = 0.5 mm,  h = 3 mm, then rotated 90° about Z
//!              → axis +X, extents X ∈ [−1.5, +1.5]
//! Both axes pass through the origin; cylinders cross at right angles.
//! ```
//!
//! The intersection of two equal perpendicular cylinders is the classical
//! **Steinmetz solid** (bicylinder).  Its exact volume is:
//!
//! ```text
//! V_∩ = (16 / 3) r³
//! ```
//!
//! With r = 0.5: V_∩ = 16/3 · 0.125 = 2/3 ≈ 0.6667 mm³
//!
//! | Operation | Formula                  | Expected (mm³) |
//! |-----------|--------------------------|----------------|
//! | A ∪ B     | 2·V_cyl − V_∩           | ≈ 3.9898       |
//! | A ∩ B     | 16r³/3                  | ≈ 0.6667       |
//! | A \ B     | V_cyl − V_∩             | ≈ 1.6616       |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_perpendicular
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use nalgebra::{Isometry3, UnitQuaternion, Vector3};

use cfd_mesh::core::scalar::{Point3r, Real};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::csg::CsgNode;
use cfd_mesh::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::io::stl;
use cfd_mesh::{IndexedMesh, analyze_normals};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (Perpendicular / Cross)");
    println!("  Union | Intersection (Steinmetz) | Difference");
    println!("  (Mesh Arrangement pipeline — axes cross at origin at 90°)");
    println!("=================================================================");

    let r: f64 = 0.5;   // both cylinder radii
    let h: f64 = 3.0;   // both cylinder heights

    let v_cyl = std::f64::consts::PI * r * r * h;

    // Steinmetz bicylinder intersection volume: 16r³/3
    let v_intersect = 16.0 / 3.0 * r * r * r;

    println!("  Cylinder A : axis +Y, base (0,{:.2},0),  r={r}, h={h}  V = {v_cyl:.4} mm³",
        -h / 2.0);
    println!("  Cylinder B : axis +X (A rotated 90° about Z), r={r}, h={h}  V = {v_cyl:.4} mm³");
    println!("  Both axes pass through origin; intersection = Steinmetz solid");
    println!("  Steinmetz V_∩ = 16r³/3 = {v_intersect:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", 2.0 * v_cyl - v_intersect);
    println!("    Intersection : {v_intersect:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_cyl - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();

    // Cylinder A: Y-axis, centred at origin.
    let cyl_a = Cylinder {
        base_center: Point3r::new(0.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
        segments: 64,
    }.build()?;

    // Cylinder B: the Cylinder primitive is Y-axis only, so we build it
    // identically to A and then rotate it into the X-axis orientation.
    // Rotating −90° about Z maps +Y → +X.
    // CsgNode::Transform applies the isometry to every vertex and normal.
    let cyl_b = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }.build()?;
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2,  // +Y → +X
        );
        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso:  Isometry3::from_parts(nalgebra::Translation3::identity(), rot),
        }.evaluate()?
    };

    println!("  Mesh built: {} + {} cylinder faces  ({} ms)",
        cyl_a.face_count(), cyl_b.face_count(), t_build.elapsed().as_millis());
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Union (A ∪ B)", &mut result, 2.0 * v_cyl - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_perpendicular_union.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_perpendicular_union.stl");
        println!();
    }

    // ── Intersection (Steinmetz solid) ────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B) — Steinmetz solid", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_perpendicular_intersection.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_perpendicular_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_cyl - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_perpendicular_difference.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_perpendicular_difference.stl");
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
