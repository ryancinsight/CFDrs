//! Export: I-WP Sphere → `outputs/primitives/iwp_sphere.stl`
//!
//! Builds an I-WP TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! and exports binary STL.
//!
//! The I-WP (Schoen 1970):
//! `F = 2(cos(kx)cos(ky) + cos(ky)cos(kz) + cos(kz)cos(kx)) − (cos(2kx)+cos(2ky)+cos(2kz)) = 0`
//!
//! Space group Im3̄m (#229), genus 4 per unit cell.  BCC topology with two
//! non-congruent interlocking labyrinths.  High genus gives excellent
//! surface-area-to-volume ratio for heat-exchanger applications.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_iwp_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{IwpSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: I-WP Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = IwpSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : I-WP (Schoen 1970 — I-graph Wrapped Package)");
    println!("  Space group: Im3̄m (#229), genus 4, BCC topology");
    println!("  Radius     : {} mm", builder.radius);
    println!(
        "  Period     : {} mm  (lattice unit-cell size)",
        builder.period
    );
    println!("  Resolution : {}  (voxels per axis)", builder.resolution);
    println!("  Building mesh...");

    let t0 = std::time::Instant::now();
    let mesh = builder.build()?;
    let elapsed = t0.elapsed();

    println!("  Vertices   : {}", mesh.vertices.len());
    println!("  Faces      : {}", mesh.faces.len());
    println!("  Elapsed    : {:?}", elapsed);

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("iwp_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
