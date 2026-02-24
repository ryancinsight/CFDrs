//! Export: Y-junction bifurcation primitives → `outputs/primitives/y_junction_*.stl`
//!
//! Demonstrates three Y-junction configurations using [`cfd_mesh::YJunction`]:
//!
//! 1. **Shallow split (30°)** — default, r=0.5 mm, 30° half-angle, typical microfluidic bifurcation
//! 2. **Wide split (45°)** — 45° half-angle, more symmetric divergence
//! 3. **Asymmetric lengths** — longer inlet, shorter branches (e.g. manifold inlet)
//!
//! All meshes are validated for watertightness (closed + oriented, χ=2)
//! and positive signed volume.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_y_junction
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::YJunction;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::application::watertight::check::check_watertight;

fn print_report(
    label: &str,
    mesh: &cfd_mesh::domain::mesh::IndexedMesh,
    stl_path: &std::path::Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let edges  = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

    println!("  [{}]", label);
    println!("    Vertices  : {}", mesh.vertices.len());
    println!("    Faces     : {}", mesh.faces.len());
    println!("    Closed    : {}", report.is_closed);
    println!("    Oriented  : {}", report.orientation_consistent);
    println!("    Watertight: {}", report.is_watertight);
    println!("    Volume    : {:.6} mm³", report.signed_volume);
    println!("    Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

    let file   = fs::File::create(stl_path)?;
    let mut wr = BufWriter::new(file);
    stl::write_binary_stl(&mut wr, &mesh.vertices, &mesh.faces)?;
    println!("    STL       : {}", stl_path.display());
    println!();
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Y-Junction  (via YJunction builder)");
    println!("=================================================================");
    println!();

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir   = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;

    // ── 1. Shallow split — 30° half-angle (default) ────────────────────────
    {
        // One inlet (r=0.5 mm, L=2 mm) splits into two 30°-diverging branches.
        // Approximate volume: π·r²·(inlet + 2·branch) = π·0.25·6 ≈ 4.71 mm³
        let yj   = YJunction::default();
        let approx_vol = PI * yj.tube_radius * yj.tube_radius
            * (yj.inlet_length + 2.0 * yj.branch_length);
        let mesh = yj.build()?;

        println!("  Approximate volume (π·r²·total_length) = {:.4} mm³", approx_vol);
        print_report(
            "Shallow split — 30° half-angle (default, r=0.5 mm, L_in=2, L_br=2)",
            &mesh,
            &out_dir.join("y_junction_30deg.stl"),
        )?;
    }

    // ── 2. Wide split — 45° half-angle ────────────────────────────────────
    {
        let yj = YJunction {
            branch_half_angle: PI / 4.0, // 45°
            tube_segments: 24,
            arm_rings: 12,
            ..YJunction::default()
        };
        let mesh = yj.build()?;
        print_report(
            "Wide split — 45° half-angle, r=0.5 mm",
            &mesh,
            &out_dir.join("y_junction_45deg.stl"),
        )?;
    }

    // ── 3. Asymmetric lengths — longer inlet, shorter branches ─────────────
    {
        // Mimics a manifold inlet feeding two short distribution channels.
        let yj = YJunction {
            tube_radius:       0.4,
            inlet_length:      4.0,
            branch_length:     1.5,
            branch_half_angle: PI / 5.0, // 36°
            tube_segments:     20,
            arm_rings:         10,
        };
        let mesh = yj.build()?;
        print_report(
            "Manifold inlet — 36° half-angle, r=0.4 mm, L_in=4.0, L_br=1.5",
            &mesh,
            &out_dir.join("y_junction_manifold.stl"),
        )?;
    }

    println!("=================================================================");
    Ok(())
}
