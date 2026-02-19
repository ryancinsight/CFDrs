//! Export: Revolution Sweep primitives → `outputs/primitives/revolution_sweep_*.stl`
//!
//! Demonstrates three profiles swept around the Y axis using
//! [`cfd_mesh::RevolutionSweep`]:
//!
//! 1. **Washer/annulus** — closed rectangular profile (r=1..2, h=0.5), full 360°
//!    Pappus volume: `R̄=1.5 × A=0.5 × 2π ≈ 4.712 mm³`
//! 2. **Quarter-turn wedge** — closed profile, swept 90°
//!    Pappus volume: `1.5 × 1 × (π/2) ≈ 2.356 mm³`
//! 3. **Half-turn elbow** — closed rectangular cross-section, swept 180°
//!    Pappus volume: `2.25 × 0.25 × π ≈ 1.767 mm³`
//!
//! All profiles are **closed** (last point == first point) so each solid is
//! fully watertight with correct Euler characteristic χ = 2.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_revolution_sweep
//! ```

use std::f64::consts::{PI, TAU};
use std::fs;
use std::io::BufWriter;

use cfd_mesh::RevolutionSweep;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Revolution Sweep  (via RevolutionSweep builder)");
    println!("=================================================================");

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;

    // ── 1. Washer — full 360° revolution of a closed rectangle ───────────
    // Closed profile (last == first): (1,0)→(2,0)→(2,0.5)→(1,0.5)→(1,0)
    // Pappus: V = 2π × R̄ × A   where R̄ = 1.5, A = 1.0 × 0.5 = 0.5
    {
        let profile = vec![
            (1.0_f64, 0.0),
            (2.0,     0.0),
            (2.0,     0.5),
            (1.0,     0.5),
            (1.0,     0.0),   // ← close the loop
        ];
        let sweep = RevolutionSweep { profile, segments: 48, angle: TAU };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

        let r_bar      = 1.5_f64;
        let area       = 1.0_f64 * 0.5_f64;
        let expected_v = TAU * r_bar * area;
        let vol        = mesh.signed_volume();
        let vol_err    = (vol - expected_v).abs() / expected_v;

        println!("  [Washer (annulus): r=1..2, h=0.5, full 360°]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (Pappus expected {:.6}, err {:.2}%)",
            vol, expected_v, vol_err * 100.0);
        println!("    Euler χ   : {:?}  (χ=0 expected for genus-1 torus topology)", report.euler_characteristic);
        let status = if report.is_watertight && vol_err < 0.02 { "PASS" } else { "FAIL" };
        println!("    Status    : {}", status);

        let stl_path = out_dir.join("revolution_sweep_washer.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    // ── 2. Quarter-turn wedge — 90° partial revolution ────────────────────
    // Closed profile: (1,0)→(2,0)→(2,1)→(1,1)→(1,0)
    // Pappus: V = (π/2) × 1.5 × 1.0 ≈ 2.356 mm³
    {
        let profile = vec![
            (1.0_f64, 0.0),
            (2.0,     0.0),
            (2.0,     1.0),
            (1.0,     1.0),
            (1.0,     0.0),   // ← close the loop
        ];
        let angle = PI / 2.0;
        let sweep = RevolutionSweep { profile, segments: 32, angle };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

        let expected_v = angle * 1.5_f64 * 1.0_f64;
        let vol        = mesh.signed_volume();
        let vol_err    = (vol - expected_v).abs() / expected_v;

        println!("  [Quarter-wedge: r=1..2, h=1, 90° sweep]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (Pappus expected {:.6}, err {:.2}%)",
            vol, expected_v, vol_err * 100.0);
        println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);
        let status = if report.is_watertight && vol_err < 0.02 { "PASS" } else { "FAIL" };
        println!("    Status    : {}", status);

        let stl_path = out_dir.join("revolution_sweep_quarter_wedge.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    // ── 3. Half-turn pipe elbow — 180° partial revolution ─────────────────
    // Closed profile (r=2..2.5, y=0..0.5): (2,0)→(2.5,0)→(2.5,0.5)→(2,0.5)→(2,0)
    // Pappus: V = π × 2.25 × 0.25 ≈ 1.767 mm³
    {
        let profile = vec![
            (2.0_f64, 0.0),
            (2.5,     0.0),
            (2.5,     0.5),
            (2.0,     0.5),
            (2.0,     0.0),   // ← close the loop
        ];
        let angle = PI;
        let sweep = RevolutionSweep { profile, segments: 32, angle };
        let mesh = sweep.build()?;

        let edges = EdgeStore::from_face_store(&mesh.faces);
        let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

        let expected_v = PI * 2.25_f64 * 0.25_f64;
        let vol        = mesh.signed_volume();
        let vol_err    = (vol - expected_v).abs() / expected_v;

        println!("  [Half-turn Elbow: r=2..2.5, h=0.5, 180° sweep]");
        println!("    Vertices  : {}", mesh.vertices.len());
        println!("    Faces     : {}", mesh.faces.len());
        println!("    Closed    : {}", report.is_closed);
        println!("    Oriented  : {}", report.orientation_consistent);
        println!("    Watertight: {}", report.is_watertight);
        println!("    Volume    : {:.6} mm³  (Pappus expected {:.6}, err {:.2}%)",
            vol, expected_v, vol_err * 100.0);
        println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);
        let status = if report.is_watertight && vol_err < 0.02 { "PASS" } else { "FAIL" };
        println!("    Status    : {}", status);

        let stl_path = out_dir.join("revolution_sweep_elbow.stl");
        let file = fs::File::create(&stl_path)?;
        let mut writer = BufWriter::new(file);
        stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;
        println!("    STL       : {}", stl_path.display());
        println!();
    }

    println!("=================================================================");
    Ok(())
}
