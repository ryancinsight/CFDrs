//! CSG Disk–Disk (Coplanar 2-D): union, intersection, and difference
//!
//! Two flat circular disks of equal radius `r = 1.0 mm`, both lying in the
//! plane `y = 0`, with centres offset by `d = r = 1.0 mm` along the X axis.
//! Because both operands are flat and coplanar the **2-D Sutherland-Hodgman
//! pipeline** (`csg::coplanar`) is used automatically — neither BSP nor the
//! 3-D Mesh-Arrangement path is invoked.
//!
//! ## Geometry
//!
//! ```text
//! Disk A : centre (−0.5, 0, 0), r = 1.0 mm   area = π mm²
//! Disk B : centre (+0.5, 0, 0), r = 1.0 mm   area = π mm²
//!   axis separation d = 1.0 mm = r  →  θ = arccos(d/2r) = arccos(0.5) = π/3
//! ```
//!
//! The intersection cross-section is a symmetric lens:
//!
//! ```text
//! θ           = π / 3
//! A_segment   = r² (θ − sin θ cos θ)
//! A_∩ = A_lens = 2 · A_segment = 2 r² (π/3 − √3/4)
//! ```
//!
//! With r = 1.0:
//!   A_segment ≈ 1.0 · (1.0472 − 0.4330) ≈ 0.6142 mm²
//!   A_∩       ≈ 2 · 0.6142              ≈ 1.2283 mm²
//!
//! | Operation | Expected (mm²)             | Pipeline   |
//! |-----------|---------------------------|------------|
//! | A ∪ B     | 2·π − A_∩ ≈ 5.0548        | Coplanar   |
//! | A ∩ B     | A_∩ ≈ 1.2283             | Coplanar   |
//! | A \ B     | π − A_∩ ≈ 1.9133         | Coplanar   |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_coplanar_disk_disk
//! ```
//!
//! STL outputs are written to `outputs/csg_coplanar/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::Point3r;
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::primitives::{Disk, PrimitiveMesh};
use cfd_mesh::io::stl;
use cfd_mesh::IndexedMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Disk–Disk (Coplanar 2-D): Union | Intersection | Difference");
    println!("  (2-D Sutherland-Hodgman pipeline — both operands flat, coplanar)");
    println!("=================================================================");

    let r: f64 = 1.0;   // both disk radii
    let d: f64 = r;     // centre separation = r  →  θ = π/3

    // Analytical lens area (circular segment × 2)
    let theta       = (d / (2.0 * r)).acos();              // π/3
    let a_seg       = r * r * (theta - theta.sin() * theta.cos());
    let a_intersect = 2.0 * a_seg;
    let a_disk      = std::f64::consts::PI * r * r;

    println!("  Disk A : centre (−0.5, 0, 0), r={r}  area = {a_disk:.4} mm²");
    println!("  Disk B : centre (+0.5, 0, 0), r={r}  area = {a_disk:.4} mm²");
    println!("  Axis sep d={d} mm (=r → θ=π/3); both in plane y=0");
    println!("  Lens area (A_∩) = {a_intersect:.4} mm²");
    println!();
    println!("  Expected areas:");
    println!("    Union        : {:.4} mm²", 2.0 * a_disk - a_intersect);
    println!("    Intersection : {a_intersect:.4} mm²");
    println!("    Difference   : {:.4} mm²", a_disk - a_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg_coplanar");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let disk_a = Disk {
        center:   Point3r::new(-d / 2.0, 0.0, 0.0),
        radius:   r,
        segments: 128,
    }.build()?;
    let disk_b = Disk {
        center:   Point3r::new( d / 2.0, 0.0, 0.0),
        radius:   r,
        segments: 128,
    }.build()?;
    println!("  Mesh built: {} + {} disk faces  ({} ms)",
        disk_a.face_count(), disk_b.face_count(), t_build.elapsed().as_millis());
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let result = csg_boolean_indexed(BooleanOp::Union, &disk_a, &disk_b)?;
        let ms = t0.elapsed().as_millis();
        report("Union (A ∪ B)", &result, 2.0 * a_disk - a_intersect, 0.02, ms);
        write_stl(&result, &out_dir.join("disk_disk_union.stl"))?;
        println!("  STL: outputs/csg_coplanar/disk_disk_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let result = csg_boolean_indexed(BooleanOp::Intersection, &disk_a, &disk_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &result, a_intersect, 0.02, ms);
        write_stl(&result, &out_dir.join("disk_disk_intersection.stl"))?;
        println!("  STL: outputs/csg_coplanar/disk_disk_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let result = csg_boolean_indexed(BooleanOp::Difference, &disk_a, &disk_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &result, a_disk - a_intersect, 0.02, ms);
        write_stl(&result, &out_dir.join("disk_disk_difference.stl"))?;
        println!("  STL: outputs/csg_coplanar/disk_disk_difference.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn mesh_area(mesh: &IndexedMesh) -> f64 {
    mesh.faces.iter().map(|f| {
        let a = mesh.vertices.position(f.vertices[0]);
        let b = mesh.vertices.position(f.vertices[1]);
        let c = mesh.vertices.position(f.vertices[2]);
        let ab = b - a;
        let ac = c - a;
        ab.cross(&ac).norm() * 0.5
    }).sum()
}

fn report(label: &str, mesh: &IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let area   = mesh_area(mesh);
    let err    = (area - expected).abs() / expected.abs().max(1e-12);
    let status = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Area       : {area:.4} mm²  (expected {expected:.4})");
    println!("    Area error : {:.2}%  [{status}]", err * 100.0);
    println!("    Elapsed    : {ms} ms");
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file  = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
