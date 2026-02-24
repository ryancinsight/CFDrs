//! CSG Sphere-Sphere: union, intersection (lens), and difference
//!
//! Demonstrates the **Mesh Arrangement CSG pipeline** for curved surfaces.
//! Two unit spheres (r = 1.0 mm) with centres 1.0 mm apart are combined
//! with all three Boolean operations.
//!
//! ## Analytical reference volumes
//!
//! For two spheres of radius r = 1 with centres separated by d = 1:
//!
//! | Operation  | Formula                          | Value (mm³) |
//! |-----------|----------------------------------|-------------|
//! | Lens (A∩B) | π/12 · (4r + d)(2r − d)²       | ≈ 1.3090    |
//! | Union (A∪B)| 2·(4π/3)r³ − V_lens              | ≈ 7.1272    |
//! | Diff (A\B) | (4π/3)r³ − V_lens               | ≈ 2.8798    |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_sphere_sphere
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::application::csg::boolean::{csg_boolean_indexed, BooleanOp};
use cfd_mesh::domain::core::scalar::{Point3r, Real};
use cfd_mesh::domain::geometry::primitives::{PrimitiveMesh, UvSphere};
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{analyze_normals, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Sphere-Sphere: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — exact curved-surface CSG)");
    println!("=================================================================");

    let r: f64 = 1.0; // sphere radius [mm]
    let d: f64 = 1.0; // centre separation [mm]

    // Analytical volumes.
    let v_sphere = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    let v_lens = std::f64::consts::PI / 12.0 * (4.0 * r + d) * (2.0 * r - d).powi(2);
    let v_union = 2.0 * v_sphere - v_lens;
    let v_diff = v_sphere - v_lens;

    println!("  r = {r:.1} mm,  d = {d:.1} mm");
    println!("  Sphere A : centre (0,0,0)   V = {v_sphere:.4} mm³");
    println!("  Sphere B : centre (1,0,0)   V = {v_sphere:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Intersection (lens) : {v_lens:.4} mm³");
    println!("    Union               : {v_union:.4} mm³");
    println!("    Difference A\\B      : {v_diff:.4} mm³");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // Build spheres with 64×32 resolution for < 5% volume error at the seam.
    let t_build = Instant::now();
    let sphere_a = UvSphere {
        radius: r,
        center: Point3r::new(0.0, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }
    .build()?;
    let sphere_b = UvSphere {
        radius: r,
        center: Point3r::new(d, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }
    .build()?;
    println!(
        "  Mesh built: {} + {} faces  ({} ms)",
        sphere_a.face_count(),
        sphere_b.face_count(),
        t_build.elapsed().as_millis()
    );
    println!();

    // ── Intersection (lens) ────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &sphere_a, &sphere_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (lens)", &mut result, v_lens, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_sphere_intersection.stl"))?;
        println!("  STL: outputs/csg/sphere_sphere_intersection.stl");
        println!();
    }

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &sphere_a, &sphere_b)?;
        let ms = t0.elapsed().as_millis();
        report("Union", &mut result, v_union, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_sphere_union.stl"))?;
        println!("  STL: outputs/csg/sphere_sphere_union.stl");
        println!();
    }

    // ── Difference A \ B ──────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &sphere_a, &sphere_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference A\\B", &mut result, v_diff, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_sphere_difference.stl"))?;
        println!("  STL: outputs/csg/sphere_sphere_difference.stl");
        println!();
    }

    println!("=================================================================");
    Ok(())
}

// ── Helpers ────────────────────────────────────────────────────────────────────

fn report(label: &str, mesh: &mut IndexedMesh, expected: f64, tol: f64, ms: u128) {
    let vol = mesh.signed_volume();
    let is_wt = mesh.is_watertight();
    let normals = analyze_normals(mesh);
    let err = (vol - expected).abs() / expected.abs().max(1e-12);
    let status = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected {expected:.4})");
    println!("    Vol error  : {:.2}%  [{status}]", err * 100.0);
    println!("    Watertight : {is_wt}");
    println!(
        "    Normals    : outward={}, inward={} ({:.1}%), degen={}",
        normals.outward_faces,
        normals.inward_faces,
        if mesh.face_count() > 0 {
            normals.inward_faces as Real / mesh.face_count() as Real * 100.0
        } else {
            0.0
        },
        normals.degenerate_faces
    );
    println!(
        "    Alignment  : mean={:.4}  min={:.4}",
        normals.face_vertex_alignment_mean, normals.face_vertex_alignment_min
    );
    println!("    Elapsed    : {ms} ms");
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}
