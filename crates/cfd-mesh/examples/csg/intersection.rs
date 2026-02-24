//! CSG Intersection: CubeA ∩ CubeB → `outputs/csg/intersection_cube_cube.stl`
//!
//! Demonstrates the Boolean intersection of two overlapping cubes.
//! CubeA is 2×2×2 mm at the origin.  CubeB is 2×2×2 mm offset by (1,1,1).
//! The two cubes share a 1×1×1 mm corner overlap region.
//! The intersection is exactly that 1×1×1 mm cube.
//!
//! Expected volume: 1 mm³  (exact, flat-face BSP)
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example csg_intersection
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
    println!("  CSG Intersection: CubeA ∩ CubeB  (shared corner = 1 mm³ cube)");
    println!("=================================================================");

    // CubeA: [0,2]³   CubeB: [1,3]³  → intersection: [1,2]³ = 1 mm³
    let expected = 1.0_f64; // 1×1×1 mm³

    println!("  Cube A : 2×2×2 mm, origin (0,0,0)");
    println!("  Cube B : 2×2×2 mm, origin (1,1,1)");
    println!("  Expected intersection: {:.4} mm³  (1×1×1 shared corner)", expected);
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    let t0 = Instant::now();

    let cube_a = Cube { origin: Point3r::new(0.0, 0.0, 0.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;
    let cube_b = Cube { origin: Point3r::new(1.0, 1.0, 1.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;

    println!("  Cube A : {} faces", cube_a.face_count());
    println!("  Cube B : {} faces", cube_b.face_count());

    let mut mesh = csg_boolean_indexed(BooleanOp::Intersection, &cube_a, &cube_b)?;

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

    // Flat-face BSP is exact — tight volume tolerance
    let vol_ok        = vol_err <= 0.02;
    let align_mean_ok = normals.face_vertex_alignment_mean >= 0.30;
    let status = if vol_ok && align_mean_ok { "PASS" } else { "FAIL" };
    println!("  Status: {} (vol_err={:.2}%, align_mean={:.4})",
        status, vol_err * 100.0, normals.face_vertex_alignment_mean);

    let stl_path = out_dir.join("intersection_cube_cube.stl");
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

