//! CSG Difference: Cube − Slot → `outputs/csg/difference_slotted_block.stl`
//!
//! A 3×3×3 mm block with a rectangular slot (channel) cut straight through it.
//!
//! - Block  (A): [0,3]³  — 27 mm³
//! - Slot   (B): [0.75, 2.25] × [0.75, 2.25] × [−1, 4]  — 1.5×1.5×5 = 11.25 mm³
//!   The slot extends ±1 mm beyond the block faces on both Z ends so it cleanly
//!   punches all the way through without a coplanar face issue.
//!
//! Overlap of B with A = 1.5 × 1.5 × 3 = 6.75 mm³  (the material removed)
//!
//! Expected volume: V_block − V_overlap = 27 − 6.75 = 20.25 mm³
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example csg_difference
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::{Point3r, Real, Vector3r};
use cfd_mesh::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::geometry::normal::triangle_normal;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::{Cube, IndexedMesh};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Difference: Block − Slot  (square channel through block)");
    println!("=================================================================");

    let v_block   = 3.0_f64.powi(3);          // 27 mm³
    let v_overlap = 1.5_f64 * 1.5 * 3.0;     // 6.75 mm³  (slot clipped to block)
    let expected  = v_block - v_overlap;       // 20.25 mm³

    println!("  Block A : 3×3×3 mm, origin (0,0,0)     V = {:.4} mm³", v_block);
    println!("  Slot  B : 1.5×1.5×5 mm, centred on XY  V_overlap = {:.4} mm³", v_overlap);
    println!("  Expected: V_block − V_overlap = {:.4} mm³", expected);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t0 = Instant::now();

    // A: 3×3×3 mm block
    let block = Cube {
        origin: Point3r::new(0.0, 0.0, 0.0),
        width: 3.0, height: 3.0, depth: 3.0,
    }.build()?;

    // B: slot that punches through from z=-1 to z=4 (1 mm beyond each face)
    //    centred on the block in XY: [0.75, 2.25]² in XY
    let slot = Cube {
        origin: Point3r::new(0.75, 0.75, -1.0),
        width: 1.5, height: 1.5, depth: 5.0,
    }.build()?;

    println!("  Block : {} faces", block.face_count());
    println!("  Slot  : {} faces", slot.face_count());

    let mut mesh = csg_boolean_indexed(BooleanOp::Difference, &block, &slot)?;

    let volume      = mesh.signed_volume();
    let is_wt       = mesh.is_watertight();
    let normals     = analyze_normals(&mesh);
    let total       = mesh.face_count();
    let inward_frac = if total > 0 { normals.inward_faces as Real / total as Real } else { 1.0 };
    let vol_err     = (volume - expected).abs() / expected.abs().max(1e-12);

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

    // 2% tolerance — cube-only BSP is numerically exact for axis-aligned faces
    let vol_ok        = vol_err <= 0.02;
    let align_mean_ok = normals.face_vertex_alignment_mean >= 0.05;
    let status = if vol_ok && align_mean_ok { "PASS" } else { "FAIL" };
    println!("  Status: {} (vol_err={:.2}%, align_mean={:.4})",
        status, vol_err * 100.0, normals.face_vertex_alignment_mean);

    let stl_path = out_dir.join("difference_slotted_block.stl");
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
