//! Export: Ring Torus → `outputs/primitives/torus.stl`
//!
//! Builds a ring torus (R=3 mm, r=1 mm) using [`cfd_mesh::Torus`],
//! validates topology, and exports binary STL.
//!
//! The torus has genus 1, so `V − E + F = 0` (Euler characteristic = 0).
//! It is still watertight (closed + consistently oriented).
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_torus
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::Torus;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Ring Torus R=3 mm, r=1 mm  (via Torus builder)");
    println!("=================================================================");

    let torus = Torus {
        major_radius: 3.0,
        minor_radius: 1.0,
        major_segments: 48,
        minor_segments: 24,
    };
    let mesh = torus.build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
    let expected_vol = 2.0 * PI * PI * 3.0 * 1.0 * 1.0;

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Closed    : {}", report.is_closed);
    println!("  Oriented  : {}", report.orientation_consistent);
    println!("  Watertight: {}", report.is_watertight);
    println!(
        "  Volume    : {:.6} mm³  (expected {:.6})",
        report.signed_volume, expected_vol
    );
    println!(
        "  Euler χ   : {:?}  (expected 0 — genus 1 torus)",
        report.euler_characteristic
    );

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("torus.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
