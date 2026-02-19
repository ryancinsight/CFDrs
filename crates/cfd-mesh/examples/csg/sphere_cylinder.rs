//! CSG Sphere–Cylinder: union, intersection, and difference
//!
//! A unit sphere (r = 1 mm) and a tall cylinder (r = 0.4 mm, h = 3 mm)
//! sharing a common Y axis.  The cylinder base is at `(0, −1.5, 0)` and the
//! top at `(0, 1.5, 0)`, extending 0.5 mm past each sphere pole.  The
//! cylinder therefore protrudes cleanly through the sphere surface at both
//! ends.  Both operands are curved, so the **Mesh Arrangement pipeline** is
//! used automatically.
//!
//! ## Geometry
//!
//! ```text
//! Sphere   A : centre (0,0,0),    r = 1 mm     V_A = 4π/3 ≈ 4.1888 mm³
//! Cylinder B : base   (0,−1.5,0), r = 0.4 mm,
//!              h = 3 mm           V_B = π·r²·h = π·0.16·3 ≈ 1.5080 mm³
//! ```
//!
//! The cylinder intersects the sphere surface at two latitude circles
//! y = ±√(1 − 0.4²) = ±√0.84 ≈ ±0.9165.  The overlap (A ∩ B) is the
//! **cylinder segment** that lies inside the sphere — a finite cylinder of
//! radius ρ and half-length h_cap:
//!
//! ```text
//! h_cap = √(R² − ρ²) = √0.84 ≈ 0.9165 mm
//! V_∩   = π ρ² · 2 h_cap  ≈  π · 0.16 · 1.8330  ≈  0.9215 mm³
//! ```
//!
//! | Operation | Expected (mm³)          | Pipeline    |
//! |-----------|------------------------|-------------|
//! | A ∪ B     | V_A + V_B − V_∩        | Arrangement |
//! | A ∩ B     | V_∩ ≈ 0.9215           | Arrangement |
//! | A \ B     | V_A − V_∩ ≈ 3.2673     | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_sphere_cylinder
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::geometry::primitives::{Cylinder, PrimitiveMesh, UvSphere};
use cfd_mesh::io::stl;
use cfd_mesh::IndexedMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Sphere–Cylinder: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — both operands curved)");
    println!("=================================================================");

    let big_r: f64  = 1.0;    // sphere radius
    let cyl_r: f64  = 0.4;    // cylinder radius — protrudes through sphere surface
    let cyl_h: f64  = 3.0;    // cylinder height — extends 0.5 mm past each pole
    let cyl_base_y  = -1.5_f64; // base Y = −1.5, top Y = 1.5

    let v_sphere = 4.0 * std::f64::consts::PI * big_r.powi(3) / 3.0;
    let v_cyl    = std::f64::consts::PI * cyl_r * cyl_r * cyl_h;

    // Analytical intersection = cylinder segment inside the sphere:
    //   h_cap = √(R²−ρ²)  — half the chord length along the axis
    //   V_∩   = π ρ² · 2 h_cap
    let h_cap = (big_r * big_r - cyl_r * cyl_r).sqrt();
    let v_intersect = std::f64::consts::PI * cyl_r * cyl_r * 2.0 * h_cap;

    println!("  Sphere   A : centre (0,0,0),    r = {big_r}  V = {v_sphere:.4} mm³");
    println!("  Cylinder B : base (0,{cyl_base_y},0), r = {cyl_r}, h = {cyl_h}  V = {v_cyl:.4} mm³");
    println!("  Overlap    : cylinder segment inside sphere   V = {v_intersect:.4} mm³");
    println!();
    println!("  Expected volumes:");
    println!("    Union        : {:.4} mm³", v_sphere + v_cyl - v_intersect);
    println!("    Intersection : {v_intersect:.4} mm³");
    println!("    Difference   : {:.4} mm³", v_sphere - v_intersect);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t_build = Instant::now();
    let sphere = UvSphere {
        radius: big_r,
        center: Point3r::new(0.0, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }.build()?;
    // Cylinder: base at (0, −1.5, 0), height 3 → top at (0, 1.5, 0).
    // Protrudes 0.5 mm past each sphere pole → surface properly intersected.
    let cylinder = Cylinder {
        base_center: Point3r::new(0.0, cyl_base_y, 0.0),
        radius: cyl_r,
        height: cyl_h,
        segments: 64,
    }.build()?;
    println!("  Mesh built: {} sphere + {} cylinder faces  ({} ms)",
        sphere.face_count(), cylinder.face_count(), t_build.elapsed().as_millis());
    println!();

    // ── Union ─────────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Union, &sphere, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report("Union (A ∪ B)", &mut result, v_sphere + v_cyl - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_cylinder_union.stl"))?;
        println!("  STL: outputs/csg/sphere_cylinder_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &sphere, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_cylinder_intersection.stl"))?;
        println!("  STL: outputs/csg/sphere_cylinder_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &sphere, &cylinder)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_sphere - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("sphere_cylinder_difference.stl"))?;
        println!("  STL: outputs/csg/sphere_cylinder_difference.stl");
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

struct NormalAnalysis {
    outward_faces: usize,
    inward_faces: usize,
    degenerate_faces: usize,
    face_vertex_alignment_mean: Real,
    face_vertex_alignment_min: Real,
}

fn analyze_normals(mesh: &IndexedMesh) -> NormalAnalysis {
    let mut centroid_sum = Vector3r::zeros();
    let mut cnt = 0usize;
    for (_, v) in mesh.vertices.iter() { centroid_sum += v.position.coords; cnt += 1; }
    let center = if cnt > 0 { Point3r::from(centroid_sum / cnt as Real) } else { Point3r::origin() };

    let mut outward = 0usize;
    let mut inward  = 0usize;
    let mut degen   = 0usize;
    let mut asum: Real = 0.0;
    let mut acnt       = 0usize;
    let mut amin: Real = 1.0;

    for face in mesh.faces.iter() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let Some(fn_) = triangle_normal(a, b, c) else { degen += 1; continue; };
        let fc  = Point3r::new((a.x+b.x+c.x)/3.0, (a.y+b.y+c.y)/3.0, (a.z+b.z+c.z)/3.0);
        let dir = fc - center;
        if dir.norm() > 1e-12 {
            if fn_.dot(&dir.normalize()) >= 0.0 { outward += 1; } else { inward += 1; }
        }
        let avg = (*mesh.vertices.normal(face.vertices[0])
            + *mesh.vertices.normal(face.vertices[1])
            + *mesh.vertices.normal(face.vertices[2])) / 3.0;
        let l = avg.norm();
        if l > 1e-12 {
            let al = fn_.dot(&(avg / l));
            asum += al; acnt += 1; amin = amin.min(al);
        }
    }
    NormalAnalysis {
        outward_faces: outward,
        inward_faces: inward,
        degenerate_faces: degen,
        face_vertex_alignment_mean: if acnt > 0 { asum / acnt as Real } else { 0.0 },
        face_vertex_alignment_min:  if acnt > 0 { amin } else { 0.0 },
    }
}
