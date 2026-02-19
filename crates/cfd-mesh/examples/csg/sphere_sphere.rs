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

use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::geometry::primitives::{PrimitiveMesh, UvSphere};
use cfd_mesh::io::stl;
use cfd_mesh::IndexedMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Sphere-Sphere: Union | Intersection | Difference");
    println!("  (Mesh Arrangement pipeline — exact curved-surface CSG)");
    println!("=================================================================");

    let r: f64 = 1.0;  // sphere radius [mm]
    let d: f64 = 1.0;  // centre separation [mm]

    // Analytical volumes.
    let v_sphere = 4.0 * std::f64::consts::PI * r.powi(3) / 3.0;
    let v_lens   = std::f64::consts::PI / 12.0 * (4.0 * r + d) * (2.0 * r - d).powi(2);
    let v_union  = 2.0 * v_sphere - v_lens;
    let v_diff   = v_sphere - v_lens;

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
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // Build spheres with 64×32 resolution for < 5% volume error at the seam.
    let t_build = Instant::now();
    let sphere_a = UvSphere {
        radius: r,
        center: Point3r::new(0.0, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }.build()?;
    let sphere_b = UvSphere {
        radius: r,
        center: Point3r::new(d, 0.0, 0.0),
        stacks: 64,
        segments: 32,
    }.build()?;
    println!("  Mesh built: {} + {} faces  ({} ms)",
        sphere_a.face_count(), sphere_b.face_count(), t_build.elapsed().as_millis());
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
    let vol     = mesh.signed_volume();
    let is_wt   = mesh.is_watertight();
    let normals = analyze_normals(mesh);
    let err     = (vol - expected).abs() / expected.abs().max(1e-12);
    let status  = if err <= tol { "PASS" } else { "FAIL" };

    println!("  ── {label} ──");
    println!("    Faces      : {}", mesh.face_count());
    println!("    Volume     : {vol:.4} mm³  (expected {expected:.4})");
    println!("    Vol error  : {:.2}%  [{status}]", err * 100.0);
    println!("    Watertight : {is_wt}");
    println!("    Normals    : outward={}, inward={} ({:.1}%), degen={}",
        normals.outward_faces, normals.inward_faces,
        if mesh.face_count() > 0 { normals.inward_faces as Real / mesh.face_count() as Real * 100.0 } else { 0.0 },
        normals.degenerate_faces);
    println!("    Alignment  : mean={:.4}  min={:.4}",
        normals.face_vertex_alignment_mean, normals.face_vertex_alignment_min);
    println!("    Elapsed    : {ms} ms");
}

fn write_stl(mesh: &IndexedMesh, path: &std::path::Path) -> Result<(), Box<dyn std::error::Error>> {
    let file   = fs::File::create(path)?;
    let mut w  = BufWriter::new(file);
    stl::write_binary_stl(&mut w, &mesh.vertices, &mesh.faces)?;
    Ok(())
}

// ── Normal analysis (shared with csg/difference.rs pattern) ─────────────────

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
    for (_, v) in mesh.vertices.iter() {
        centroid_sum += v.position.coords;
        cnt += 1;
    }
    let mesh_center = if cnt > 0 {
        Point3r::from(centroid_sum / cnt as Real)
    } else {
        Point3r::origin()
    };

    let mut outward_faces   = 0usize;
    let mut inward_faces    = 0usize;
    let mut degenerate_faces = 0usize;
    let mut align_sum: Real  = 0.0;
    let mut align_count      = 0usize;
    let mut align_min: Real  = 1.0;

    for face in mesh.faces.iter() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let Some(face_n) = triangle_normal(a, b, c) else {
            degenerate_faces += 1;
            continue;
        };
        let fc = Point3r::new((a.x + b.x + c.x) / 3.0,
                               (a.y + b.y + c.y) / 3.0,
                               (a.z + b.z + c.z) / 3.0);
        let to_face = fc - mesh_center;
        if to_face.norm() > 1e-12 {
            if face_n.dot(&to_face.normalize()) >= 0.0 {
                outward_faces += 1;
            } else {
                inward_faces += 1;
            }
        }
        let avg_n = (*mesh.vertices.normal(face.vertices[0])
            + *mesh.vertices.normal(face.vertices[1])
            + *mesh.vertices.normal(face.vertices[2])) / 3.0;
        let l = avg_n.norm();
        if l > 1e-12 {
            let al = face_n.dot(&(avg_n / l));
            align_sum   += al;
            align_count += 1;
            align_min    = align_min.min(al);
        }
    }

    NormalAnalysis {
        outward_faces,
        inward_faces,
        degenerate_faces,
        face_vertex_alignment_mean: if align_count > 0 { align_sum / align_count as Real } else { 0.0 },
        face_vertex_alignment_min:  if align_count > 0 { align_min } else { 0.0 },
    }
}
