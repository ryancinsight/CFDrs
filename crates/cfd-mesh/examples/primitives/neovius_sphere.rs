//! Export: Neovius Sphere → `outputs/primitives/neovius_sphere.stl`
//!
//! Builds a Neovius TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! and exports binary STL.
//!
//! The Neovius surface: `F = 3(cos(kx)+cos(ky)+cos(kz)) + 4·cos(kx)·cos(ky)·cos(kz) = 0`
//!
//! Space group Im3̄m (#229), genus 9 per unit cell.  Bicontinuous labyrinthine
//! network with two congruent interlocking channel systems.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_neovius_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{NeoviusSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Neovius Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = NeoviusSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : F = 3(cos(kx)+cos(ky)+cos(kz)) + 4·cos(kx)cos(ky)cos(kz)");
    println!("  Space group: Im3̄m (#229), genus 9 per cell");
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
    let stl_path = out_dir.join("neovius_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
