//! CSG Cube–Sphere: union, intersection, and difference
//!
//! A 2×2×2 mm axis-aligned cube and a unit sphere (r = 1 mm) whose centre
//! sits 0.5 mm inside the cube's +X face at `(1.5, 1.0, 1.0)`.
//! The sphere protrudes 0.5 mm past the face — a spherical cap sticks out of
//! the +X side.  Because one operand is curved the **Mesh Arrangement
//! pipeline** is used automatically.
//!
//! ## Geometry
//!
//! ```text
//! Cube   A : [0, 2]³ mm                   V_A = 8 mm³
//! Sphere B : centre (1.5, 1, 1), r = 1    V_B = 4π/3 ≈ 4.1888 mm³
//! ```
//!
//! The +X face (x = 2) intersects the sphere, cutting off a spherical cap of
//! height h = R − d = 1 − 0.5 = 0.5 mm.  The cap that sticks **outside** the
//! cube has volume:
//!
//! ```text
//! V_cap = (π h²/3)(3R − h) = (π · 0.25 / 3)(2.5) = 5π/24 ≈ 0.6545 mm³
//! ```
//!
//! The overlap (inside the cube) is:
//!
//! ```text
//! V_∩ = V_B − V_cap ≈ 3.5343 mm³
//! ```
//!
//! | Operation | Expected (mm³)          | Pipeline    |
//! |-----------|------------------------|-------------|
//! | A ∪ B     | V_A + V_B − V_∩        | Arrangement |
//! | A ∩ B     | V_∩ ≈ 3.5343           | Arrangement |
//! | A \ B     | V_A − V_∩ ≈ 4.4657     | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cube_sphere
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::primitives::{PrimitiveMesh, UvSphere};
use cfd_mesh::io::stl;
use cfd_mesh::{Cube, IndexedMesh, analyze_normals};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cube–Sphere: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — curved surface)");
    println!("=================================================================");

    let r: f64 = 1.0;
    let d: f64 = 0.5; // distance from sphere centre to +X face (x = 2 − 1.5 = 0.5)
    let v_cube   = 8.0_f64;
    let v_sphere = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    // Spherical cap outside cube: h = R − d = 0.5
    let h_cap  = r - d;
    let v_cap  = std::f64::consts::PI * h_cap * h_cap / 3.0 * (3.0 * r - h_cap);
    let v_overlap = v_sphere - v_cap;  // portion inside the cube

    println!("  Cube   A : [0,2]³ mm                         V = {v_cube:.4} mm³");
    println!("  Sphere B : centre (1.5,1,1), r = {r:.1} mm   V = {v_sphere:.4} mm³");
    println!("  Cap outside cube (h={h_cap})                  V = {v_cap:.4} mm³");
    println!("  Overlap (sphere inside cube)                   V = {v_overlap:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", v_cube + v_sphere - v_overlap);
    println!("    Intersection : {v_overlap:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_cube - v_overlap);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let cube = Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 2.0, height: 2.0, depth: 2.0,
    }.build()?;
    let sphere = UvSphere {
        radius: r,
        center: Point3r::new(1.5, 1.0, 1.0),
        stacks: 64,
        segments: 32,
    }.build()?;
    println!("  Mesh built: {} cube + {} sphere faces  ({} ms)",
        cube.face_count(), sphere.face_count(), t_build.elapsed().as_millis());
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &cube, &sphere)?;
        let ms = t0.elapsed().as_millis();
        report("Union (A ∪ B)", &mut result, v_cube + v_sphere - v_overlap, 0.05, ms);
        write_stl(&result, &out_dir.join("cube_sphere_union.stl"))?;
        println!("  STL: outputs/csg/cube_sphere_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cube, &sphere)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_overlap, 0.05, ms);
        write_stl(&result, &out_dir.join("cube_sphere_intersection.stl"))?;
        println!("  STL: outputs/csg/cube_sphere_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cube, &sphere)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_cube - v_overlap, 0.05, ms);
        write_stl(&result, &out_dir.join("cube_sphere_difference.stl"))?;
        println!("  STL: outputs/csg/cube_sphere_difference.stl");
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

