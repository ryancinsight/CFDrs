//! CSG Cylinder–Cylinder (T-Junction): union, intersection, and difference
//!
//! Two **equal-radius** cylinders arranged as a **T**: the stem runs along +Y
//! and the crossbar runs along +X, with the crossbar centred at the **top end**
//! of the stem so that their axes meet at a right angle at one endpoint.
//!
//! ## Geometry
//!
//! ```text
//!          ┌──────────────────────────┐
//!          │   crossbar  (axis +X)    │   y = h/2
//!          └──────────────────────────┘
//!                      │
//!                      │  stem (axis +Y)
//!                      │
//!                      └  y = -h/2
//!
//! Stem    (A): base (0, −h/2, 0), axis +Y, r = 0.5, h = 3 mm
//! Crossbar(B): axis +X, centred at (0, +h/2, 0), r = 0.5, h = 3 mm
//!              → built as Y-axis cylinder, translated to top of stem,
//!                then rotated −90° about Z (+Y → +X)
//! ```
//!
//! ## Volume analysis
//!
//! The crossbar sits at the very top of the stem.  Their overlap is a
//! half-Steinmetz solid: the portion of each cylinder inside the other,
//! measured at the junction.  Because the crossbar axis passes through the
//! **rim** of the stem cap (not its centre), the intersection is exactly
//! **half** the full Steinmetz bicylinder:
//!
//! ```text
//! V_full_Steinmetz = 16r³/3
//! V_∩              = (16r³/3) / 2  =  8r³/3
//! ```
//!
//! With r = 0.5:
//!   V_∩ = 8 · 0.125 / 3 = 1/3 ≈ 0.3333 mm³
//!
//! | Operation | Formula            | Expected (mm³) |
//! |-----------|--------------------|----------------|
//! | A ∪ B     | V_A + V_B − V_∩   | ≈ 4.3791       |
//! | A ∩ B     | 8r³/3              | ≈ 0.3333       |
//! | A \ B     | V_A − V_∩          | ≈ 2.0229       |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_t_junction
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use nalgebra::{Isometry3, Translation3, UnitQuaternion, Vector3};

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::application::csg::CsgNode;
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (T-Junction)");
    println!("  Union | Intersection | Difference");
    println!("  (crossbar axis meets stem top-cap at 90°)");
    println!("=================================================================");

    let r: f64 = 0.5; // both cylinder radii
    let h: f64 = 3.0; // both cylinder heights

    let v_cyl = std::f64::consts::PI * r * r * h;

    // Crossbar axis passes through the top rim of the stem (y = +h/2).
    // The overlap is half the Steinmetz bicylinder: V_∩ = 8r³/3.
    let v_intersect = 8.0 / 3.0 * r * r * r;

    println!(
        "  Stem    (A): base (0,{:.2},0), axis +Y, r={r}, h={h}  V = {v_cyl:.4} mm³",
        -h / 2.0
    );
    println!(
        "  Crossbar(B): axis +X, centred at (0,{:.2},0), r={r}, h={h}  V = {v_cyl:.4} mm³",
        h / 2.0
    );
    println!("  Crossbar meets stem at its top cap → half-Steinmetz overlap");
    println!("  V_∩ = 8r³/3 = {v_intersect:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", 2.0 * v_cyl - v_intersect);
    println!("    Intersection : {v_intersect:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_cyl - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();

    // Stem (A): Y-axis, centred at origin.
    let cyl_a = Cylinder {
        base_center: Point3r::new(0.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
        segments: 64,
    }
    .build()?;

    // Crossbar (B): Y-axis cylinder translated so its centre sits at the top
    // of the stem (0, +h/2, 0), then rotated −90° about Z (+Y → +X).
    // The Isometry3 = rotation ∘ translation (translation applied first in
    // body frame, then rotated into world frame).
    let cyl_b = {
        let raw = Cylinder {
            base_center: Point3r::new(0.0, -h / 2.0, 0.0),
            radius: r,
            height: h,
            segments: 64,
        }
        .build()?;

        // Rotate −90° about Z: +Y → +X (crossbar axis becomes +X).
        let rot = UnitQuaternion::<Real>::from_axis_angle(
            &Vector3::z_axis(),
            -std::f64::consts::FRAC_PI_2, // +Y → +X
        );
        // After rotation the cylinder is already centred at origin on the X axis.
        // Translate upward by +h/2 so the crossbar sits at the top of the stem.
        let translation = Translation3::new(0.0, h / 2.0, 0.0);
        let iso = Isometry3::from_parts(translation, rot);

        CsgNode::Transform {
            node: Box::new(CsgNode::Leaf(raw)),
            iso,
        }
        .evaluate()?
    };

    println!(
        "  Mesh built: {} + {} cylinder faces  ({} ms)",
        cyl_a.face_count(),
        cyl_b.face_count(),
        t_build.elapsed().as_millis()
    );
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Union (A ∪ B)",
            &mut result,
            2.0 * v_cyl - v_intersect,
            0.05,
            ms,
        );
        write_stl(
            &result,
            &out_dir.join("cylinder_cylinder_t_junction_union.stl"),
        )?;
        println!("  STL: outputs/csg/cylinder_cylinder_t_junction_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(
            &result,
            &out_dir.join("cylinder_cylinder_t_junction_intersection.stl"),
        )?;
        println!("  STL: outputs/csg/cylinder_cylinder_t_junction_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report(
            "Difference (A \\ B)",
            &mut result,
            v_cyl - v_intersect,
            0.05,
            ms,
        );
        write_stl(
            &result,
            &out_dir.join("cylinder_cylinder_t_junction_difference.stl"),
        )?;
        println!("  STL: outputs/csg/cylinder_cylinder_t_junction_difference.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let is_wt = mesh.is_watertight();
    let n = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let status = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected {expected:.4})");
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
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
