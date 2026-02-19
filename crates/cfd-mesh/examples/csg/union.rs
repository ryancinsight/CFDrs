//! CSG Union: CubeA ∪ CubeB → `outputs/csg/union_cube_cube.stl`
//!
//! Demonstrates the Boolean union of two overlapping cubes.
//! CubeA is 2×2×2 mm at the origin.  CubeB is 2×2×2 mm offset by (1,0,0) so
//! the two cubes share a 1×2×2 = 4 mm² overlap slab.  The result is an
//! L-shaped rectangular solid with volume 8 + 8 − 4 = 12 mm³.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example csg_union
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::Cube;
use cfd_mesh::IndexedMesh;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Union: CubeA ∪ CubeB  (overlapping slabs → L-bar)");
    println!("=================================================================");

    // CubeA: [0,2]×[0,2]×[0,2]   CubeB: [1,3]×[0,2]×[0,2]
    // Overlap: [1,2]×[0,2]×[0,2] = 1×2×2 = 4 mm³
    let v_a      = 2.0_f64.powi(3); // 8 mm³
    let v_b      = 2.0_f64.powi(3); // 8 mm³
    let v_overlap = 1.0 * 2.0 * 2.0; // 4 mm³
    let expected  = v_a + v_b - v_overlap; // 12 mm³

    println!("  Cube A : 2×2×2 mm, origin (0,0,0)  V = {:.4} mm³", v_a);
    println!("  Cube B : 2×2×2 mm, origin (1,0,0)  V = {:.4} mm³  overlap = {:.4} mm³", v_b, v_overlap);
    println!("  Expected union: {:.4} mm³", expected);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t0 = Instant::now();

    let cube_a = Cube { origin: Point3r::new(0.0, 0.0, 0.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;
    let cube_b = Cube { origin: Point3r::new(1.0, 0.0, 0.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;

    println!("  Cube A : {} faces", cube_a.face_count());
    println!("  Cube B : {} faces", cube_b.face_count());

    let mut mesh = csg_boolean_indexed(BooleanOp::Union, &cube_a, &cube_b)?;

    let volume    = mesh.signed_volume();
    let is_wt     = mesh.is_watertight();
    let normals   = analyze_normals(&mesh);
    let total     = mesh.face_count();
    let inward_frac = if total > 0 { normals.inward_faces as Real / total as Real } else { 1.0 };
    let vol_err   = (volume - expected).abs() / expected.abs().max(1e-12);

    println!("  Result : {} faces", total);
    println!();
    println!("  Volume        : {:.4} mm³  (expected {:.4})", volume, expected);
    println!("  Volume error  : {:.2}%", vol_err * 100.0);
    println!("  Watertight    : {}", is_wt);
    println!("  Normal analysis:");
    println!("    outward={}, inward={} ({:.1}%), degenerate={}",
        normals.outward_faces, normals.inward_faces, inward_frac * 100.0,
        normals.degenerate_faces);
    println!("    face↔vertex alignment: mean={:.4}, min={:.4}",
        normals.face_vertex_alignment_mean, normals.face_vertex_alignment_min);

    // Flat-face BSP is exact — tight volume tolerance; all faces outward
    let vol_ok        = vol_err <= 0.02;
    let align_mean_ok = normals.face_vertex_alignment_mean >= 0.30;
    let status = if vol_ok && align_mean_ok { "PASS" } else { "FAIL" };
    println!("  Status: {} (vol_err={:.2}%, align_mean={:.4})",
        status, vol_err * 100.0, normals.face_vertex_alignment_mean);

    let stl_path = out_dir.join("union_cube_cube.stl");
    {
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
    }
    println!("  STL     : {}", stl_path.display());
    println!("  Elapsed : {} ms", t0.elapsed().as_millis());
    println!("=================================================================");
    Ok(())
}

// =============================================================================
// Normal analysis
// =============================================================================

#[derive(Debug)]
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
    let mesh_center = if cnt > 0 { Point3r::from(centroid_sum / cnt as Real) } else { Point3r::origin() };

    let mut outward_faces = 0usize;
    let mut inward_faces = 0usize;
    let mut degenerate_faces = 0usize;
    let mut align_sum: Real = 0.0;
    let mut align_count = 0usize;
    let mut align_min: Real = 1.0;

    for face in mesh.faces.iter() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        let Some(face_n) = triangle_normal(a, b, c) else { degenerate_faces += 1; continue; };
        let fc = Point3r::new((a.x+b.x+c.x)/3.0, (a.y+b.y+c.y)/3.0, (a.z+b.z+c.z)/3.0);
        let to_face = fc - mesh_center;
        if to_face.norm() > 1e-12 {
            if face_n.dot(&to_face.normalize()) >= 0.0 { outward_faces += 1; } else { inward_faces += 1; }
        }
        let avg_n = (*mesh.vertices.normal(face.vertices[0])
            + *mesh.vertices.normal(face.vertices[1])
            + *mesh.vertices.normal(face.vertices[2])) / 3.0;
        let l = avg_n.norm();
        if l > 1e-12 {
            let al = face_n.dot(&(avg_n / l));
            align_sum += al; align_count += 1; align_min = align_min.min(al);
        }
    }
    NormalAnalysis {
        outward_faces, inward_faces, degenerate_faces,
        face_vertex_alignment_mean: if align_count > 0 { align_sum / align_count as Real } else { 0.0 },
        face_vertex_alignment_min:  if align_count > 0 { align_min } else { 0.0 },
    }
}
