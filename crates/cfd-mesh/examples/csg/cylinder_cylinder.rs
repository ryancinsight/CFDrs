//! CSG Cylinder–Cylinder: union, intersection, and difference
//!
//! Two Y-axis cylinders of equal radius (r = 0.6 mm) with parallel but
//! offset axes (Δx = 0.6 mm = r).  Cylinder A is h = 3 mm tall; Cylinder B
//! is taller (h = 4 mm) and centred so it extends 0.5 mm past both caps of A.
//! The height mismatch avoids coplanar cap faces.  Both operands are curved,
//! so the **Mesh Arrangement pipeline** is used automatically.
//!
//! ## Geometry
//!
//! ```text
//! Cylinder A : base (−0.3, −1.5, 0), r = 0.6, h = 3  →  Y ∈ [−1.5, 1.5]
//! Cylinder B : base (+0.3, −2.0, 0), r = 0.6, h = 4  →  Y ∈ [−2.0, 2.0]
//!   axis separation d = 0.6 mm = r  →  θ = arccos(d / 2r) = arccos(0.5) = π/3
//! ```
//!
//! Because B fully spans A in Y, the intersection is the cylindrical lens
//! cross-section extruded over A's full height:
//!
//! ```text
//! θ         = π / 3
//! A_segment = r² (θ − sin θ cos θ)  =  r² (π/3 − √3/4)
//! A_lens    = 2 · A_segment
//! V_∩       = A_lens · h_A  =  2 · h_A · r² · (π/3 − √3/4)
//! ```
//!
//! With r = 0.6, h_A = 3:
//!   A_lens ≈ 2 · 0.36 · (1.0472 − 0.4330) ≈ 0.4422 mm²
//!   V_∩    ≈ 0.4422 · 3 ≈ 1.3267 mm³
//!
//! | Operation | Expected (mm³)            | Pipeline    |
//! |-----------|--------------------------|-------------|
//! | A ∪ B     | V_A + V_B − V_∩          | Arrangement |
//! | A ∩ B     | V_∩ ≈ 1.3267             | Arrangement |
//! | A \ B     | V_A − V_∩ ≈ 2.0662       | Arrangement |
//!
//! ## Run
//!
//! ```sh
//! cargo run -p cfd-mesh --example csg_cylinder_cylinder
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
    println!("  CSG Cylinder–Cylinder: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — both operands curved)");
    println!("=================================================================");

    let r: f64  = 0.6;   // both cylinder radii
    let h_a: f64 = 3.0;  // Cylinder A height
    let h_b: f64 = 4.0;  // Cylinder B height (taller → B spans A in Y, no coplanar caps)
    let d: f64  = r;     // axis separation = r  →  θ = π/3

    let v_a = std::f64::consts::PI * r * r * h_a;
    let v_b = std::f64::consts::PI * r * r * h_b;

    // Cylindrical-lens cross-section (d = r):
    //   θ = arccos(d / 2r) = arccos(0.5) = π/3
    //   A_segment = r²(θ − sin θ · cos θ)
    //   V_∩ = 2 · h_A · A_segment   (B spans A's full height)
    let theta = (d / (2.0 * r)).acos();                      // π/3
    let a_seg = r * r * (theta - theta.sin() * theta.cos());
    let v_intersect = 2.0 * h_a * a_seg;

    println!("  Cylinder A : base (−0.3,−1.5,0), r={r}, h={h_a}  V = {v_a:.4} mm³");
    println!("  Cylinder B : base (+0.3,−2.0,0), r={r}, h={h_b}  V = {v_b:.4} mm³");
    println!("  Axis separation d={d} mm (= r → θ=π/3)");
    println!("  Overlap    : cylindrical lens over A's height      V = {v_intersect:.4} mm³");
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
    // A: Y ∈ [−1.5, 1.5]
    let cyl_a = Cylinder {
        base_center: Point3r::new(-d / 2.0, -h_a / 2.0, 0.0),
        radius: r,
        height: h_a,
        segments: 64,
    }.build()?;
    // B: Y ∈ [−2.0, 2.0] — extends 0.5 mm past each cap of A
    let cyl_b = Cylinder {
        base_center: Point3r::new( d / 2.0, -h_b / 2.0, 0.0),
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
        write_stl(&result, &out_dir.join("cylinder_cylinder_union.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_union.stl");
        println!();
    }

    // ── Intersection ──────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Intersection, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Intersection (A ∩ B)", &mut result, v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_intersection.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_intersection.stl");
        println!();
    }

    // ── Difference ────────────────────────────────────────────────────────────
    {
        let t0 = Instant::now();
        let mut result = csg_boolean_indexed(BooleanOp::Difference, &cyl_a, &cyl_b)?;
        let ms = t0.elapsed().as_millis();
        report("Difference (A \\ B)", &mut result, v_a - v_intersect, 0.05, ms);
        write_stl(&result, &out_dir.join("cylinder_cylinder_difference.stl"))?;
        println!("  STL: outputs/csg/cylinder_cylinder_difference.stl");
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
