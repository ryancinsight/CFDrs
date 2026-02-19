//! CSG Cylinder–Cylinder (Symmetric): union, intersection, and difference
//!
//! Two **equal-radius, equal-height** Y-axis cylinders (r = 0.6 mm, h = 3 mm)
//! whose axes are parallel but offset by exactly one radius (Δx = 0.6 mm = r).
//! The cylinders partially overlap along their full height, forming a
//! **cylindrical lens** cross-section.  Both operands are curved and their cap
//! faces are **exactly coplanar** — this exercises the coplanar-face
//! classification path in the Mesh Arrangement pipeline.
//!
//! The symmetric case (equal heights) is the harder scenario because all four
//! cap planes (A-top, B-top, A-bottom, B-bottom) are shared at y = ±1.5.
//! See `cylinder_cylinder_asymmetric.rs` for the case where caps are offset.
//!
//! ## Geometry
//!
//! ```text
//! Cylinder A : base (−0.3, −1.5, 0), r = 0.6, h = 3 mm  →  Y ∈ [−1.5, 1.5]
//! Cylinder B : base (+0.3, −1.5, 0), r = 0.6, h = 3 mm  →  Y ∈ [−1.5, 1.5]
//!   axis separation d = 0.6 mm = r  →  θ = arccos(d / 2r) = arccos(0.5) = π/3
//!   cap planes exactly coplanar at y = ±1.5
//! ```
//!
//! The intersection cross-section is a symmetric lens of two equal circular
//! segments.  For d = r:
//!
//! ```text
//! θ         = π / 3
//! A_segment = r² (θ − sin θ cos θ)  =  r² (π/3 − √3/4)
//! A_lens    = 2 · A_segment
//! V_∩       = A_lens · h  =  2 h r² (π/3 − √3/4)
//! ```
//!
//! With r = 0.6, h = 3:
//!   A_lens ≈ 2 · 0.36 · (1.0472 − 0.4330) ≈ 0.4422 mm²
//!   V_∩    ≈ 0.4422 · 3 ≈ 1.3267 mm³
//!
//! | Operation | Expected (mm³)          | Pipeline    |
//! |-----------|------------------------|-------------|
//! | A ∪ B     | 2·V_cyl − V_∩          | Arrangement |
//! | A ∩ B     | V_∩ ≈ 1.3267           | Arrangement |
//! | A \ B     | V_cyl − V_∩ ≈ 2.0662   | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder_symmetric
//! ```
//!
//! STL outputs are written to `outputs/csg/`.

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::geometry::primitives::{Cylinder, PrimitiveMesh};
use cfd_mesh::io::stl;
use cfd_mesh::IndexedMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Cylinder–Cylinder (Symmetric): Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — equal heights, exactly coplanar caps)");
    println!("=================================================================");

    let r: f64 = 0.6;   // both cylinder radii
    let h: f64 = 3.0;   // both cylinder heights (equal → coplanar cap planes)
    let d: f64 = r;     // axis separation = r  →  θ = π/3

    let v_cyl = std::f64::consts::PI * r * r * h;

    // Cylindrical-lens cross-section (d = r  →  θ = π/3):
    //   A_segment = r²(θ − sin θ · cos θ)
    //   V_∩ = 2 · h · A_segment
    let theta  = (d / (2.0 * r)).acos();              // π/3
    let a_seg  = r * r * (theta - theta.sin() * theta.cos());
    let v_intersect = 2.0 * h * a_seg;

    println!("  Cylinder A : base (−0.3,−1.5,0), r={r}, h={h}  V = {v_cyl:.4} mm³");
    println!("  Cylinder B : base (+0.3,−1.5,0), r={r}, h={h}  V = {v_cyl:.4} mm³");
    println!("  Axis sep d={d} mm (=r → θ=π/3); caps coplanar at y=±1.5");
    println!("  Overlap    : cylindrical lens cross-section      V = {v_intersect:.4} mm³");
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
    let cyl_a = Cylinder {
        base_center: Point3r::new(-d / 2.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
        segments: 64,
    }.build()?;
    let cyl_b = Cylinder {
        base_center: Point3r::new( d / 2.0, -h / 2.0, 0.0),
        radius: r,
        height: h,
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
        report("Union (A ∪ B)", &mut result, 2.0 * v_cyl - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_symmetric_union.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_symmetric_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_symmetric_intersection.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_symmetric_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_cyl - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_symmetric_difference.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_symmetric_difference.stl");
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
