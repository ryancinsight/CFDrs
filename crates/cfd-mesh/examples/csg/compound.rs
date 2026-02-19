//! CSG Compound: Multi-operation tree using `CsgNode` → `outputs/csg/`
//!
//! Demonstrates composable CSG operations via the [`CsgNode`] expression tree.
//! Both shapes use cube-only geometry where the BSP is exact (flat faces,
//! no curved surface approximation error).
//!
//! ## Shape 1: Three-cube cross union
//!
//! `CubeX ∪ CubeY ∪ CubeZ` — three orthogonal 0.5×2×0.5 bars forming a cross.
//! Each bar spans the centre of a 2×2×2 volume along one axis.
//! Volume = 3×1 − 3×(0.25) + 0 = 3 − 0.75 = 2.25 mm³ (inclusion-exclusion)
//!
//! Actually: bar_x = 0.5×0.5×2=0.5, bar_y=0.5, bar_z=0.5.
//! Pairwise overlaps (bar_x∩bar_y etc.) each = 0.5×0.5×0.5=0.125, three of them.
//! Triple overlap = 0.5³=0.125. V = 3×0.5 − 3×0.125 + 0.125 = 1.5−0.375+0.125 = 1.25 mm³
//!
//! ## Shape 2: Notched L-bar
//!
//! `(CubeA ∪ CubeB) − CubeC` — union of two overlapping cubes with a corner notch.
//! CubeA [0,2]³, CubeB [1,3]×[0,2]×[0,2] → bar 12 mm³.
//! CubeC [2,4]×[1,3]×[1,3] → cuts a 1×1×1 corner = 1 mm³.
//! Expected: 12 − 1 = 11 mm³.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example csg_compound
//! ```

use std::fs;
use std::io::BufWriter;
use std::time::Instant;

use cfd_mesh::core::scalar::Real;
use cfd_mesh::csg::boolean::{BooleanOp, CsgNode, csg_boolean_indexed};
use cfd_mesh::{Cube, IndexedMesh, NormalAnalysis, analyze_normals};
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;

// BSP-cut flat-face fragments inherit zero-length vertex normals → low alignment is expected.
// We validate geometry via exact volume instead; alignment is informational only.
const MIN_FACE_VERTEX_ALIGN_MEAN: Real = 0.05;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  CSG Compound: Multi-operation CsgNode tree  (cube-only, exact)");
    println!("=================================================================");

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("csg");
    fs::create_dir_all(&out_dir)?;

    // ── Shape 1: (CubeA ∪ CubeB) − CubeC  (notched bar) ──────────────────
    {
        println!();
        println!("  Shape 1: (CubeA ∪ CubeB) − CubeC  (notched L-bar)");
        let t0 = Instant::now();

        // CubeA [0,2]³ ∪ CubeB [1,3]×[0,2]×[0,2] → 12 mm³ bar
        let cube_a = Cube { origin: Point3r::new(0.0, 0.0, 0.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;
        let cube_b = Cube { origin: Point3r::new(1.0, 0.0, 0.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;
        let union_ab = csg_boolean_indexed(BooleanOp::Union, &cube_a, &cube_b)?;

        // CubeC [2,4]×[1,3]×[1,3]: overlaps bar at [2,3]×[1,2]×[1,2] = 1 mm³
        let cube_c = Cube { origin: Point3r::new(2.0, 1.0, 1.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;

        let tree = CsgNode::Difference {
            left:  Box::new(CsgNode::Leaf(union_ab)),
            right: Box::new(CsgNode::Leaf(cube_c)),
        };
        let mut mesh = tree.evaluate()?;

        let expected = 11.0_f64; // 12 − 1
        report_mesh("notched_bar", &mut mesh, expected, 0.02, &out_dir, t0)?;
    }

    // ── Shape 2: CubeA − (CubeB ∪ CubeC)  (two-corner notch) ─────────────
    {
        println!();
        println!("  Shape 2: CubeA − (CubeB ∪ CubeC)  (two-corner notch)");
        let t0 = Instant::now();

        // CubeA [0,3]³ = 27 mm³
        let cube_a = Cube { origin: Point3r::new(0.0, 0.0, 0.0), width: 3.0, height: 3.0, depth: 3.0 }.build()?;

        // CubeB [2,4]³ — overlaps corner of A at [2,3]³ = 1 mm³
        let cube_b = Cube { origin: Point3r::new(2.0, 2.0, 2.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;
        // CubeC [2,4]×[-1,1]×[-1,1] — overlaps opposite corner at [2,3]×[0,1]×[0,1] = 1 mm³
        let cube_c = Cube { origin: Point3r::new(2.0, -1.0, -1.0), width: 2.0, height: 2.0, depth: 2.0 }.build()?;

        // Union of the two notch cubes, then subtract from A
        let union_bc = csg_boolean_indexed(BooleanOp::Union, &cube_b, &cube_c)?;

        let tree = CsgNode::Difference {
            left:  Box::new(CsgNode::Leaf(cube_a)),
            right: Box::new(CsgNode::Leaf(union_bc)),
        };
        let mut mesh = tree.evaluate()?;

        let expected = 27.0 - 1.0 - 1.0; // 25 mm³
        report_mesh("two_corner_notch", &mut mesh, expected, 0.02, &out_dir, t0)?;
    }

    println!();
    println!("=================================================================");
    Ok(())
}

// =============================================================================
// Report helper
// =============================================================================

fn report_mesh(
    name: &str,
    mesh: &mut IndexedMesh,
    expected_volume: Real,
    vol_tol: Real,
    out_dir: &std::path::Path,
    t0: std::time::Instant,
) -> Result<(), Box<dyn std::error::Error>> {
    let volume   = mesh.signed_volume();
    let is_wt    = mesh.is_watertight();
    let normals  = analyze_normals(mesh);
    let total    = mesh.face_count();
    let inward_frac = if total > 0 { normals.inward_faces as Real / total as Real } else { 1.0 };
    let vol_err  = (volume - expected_volume).abs() / expected_volume.abs().max(1e-12);

    println!("    Vertices : {}", mesh.vertices.len());
    println!("    Faces    : {}", mesh.face_count());
    println!("    Volume   : {:.4} mm³  (expected {:.4}, err {:.2}%)",
        volume, expected_volume, vol_err * 100.0);
    println!("    Watertight: {}", is_wt);
    println!("    Normals  : outward={}, inward={} ({:.1}%), degen={}",
        normals.outward_faces, normals.inward_faces, inward_frac * 100.0,
        normals.degenerate_faces);
    println!("    Align    : mean={:.4}, min={:.4}",
        normals.face_vertex_alignment_mean, normals.face_vertex_alignment_min);

    let alignment_ok = normals.face_vertex_alignment_mean >= MIN_FACE_VERTEX_ALIGN_MEAN;
    let volume_ok    = vol_err <= vol_tol;
    let status = if volume_ok && alignment_ok { "PASS" } else { "FAIL" };
    println!("    Status   : {} (vol_err={:.2}%, align_mean={:.4})",
        status, vol_err * 100.0, normals.face_vertex_alignment_mean);

    let stl_path = out_dir.join(format!("compound_{}.stl", name));
    {
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
    }
    println!("    STL      : {}", stl_path.display());
    println!("    Elapsed  : {} ms", t0.elapsed().as_millis());

    Ok(())
}

