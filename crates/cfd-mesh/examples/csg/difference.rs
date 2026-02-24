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

use cfd_mesh::domain::core::scalar::Real;
use cfd_mesh::application::csg::boolean::{BooleanOp, csg_boolean_indexed};
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::{Cube, IndexedMesh, NormalAnalysis, analyze_normals};

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

    // 2% volume tolerance.
    // Note: face-to-centroid alignment is not a valid correctness metric for
    // non-convex solids (e.g., slotted blocks where slot-wall normals point
    // inward relative to the overall centroid but are geometrically correct).
    // Watertightness + volume accuracy are the appropriate checks here.
    let vol_ok = vol_err <= 0.02;
    let status = if vol_ok && is_wt { "PASS" } else { "FAIL" };
    println!("  Status: {} (vol_err={:.2}%, watertight={})",
        status, vol_err * 100.0, is_wt);

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

