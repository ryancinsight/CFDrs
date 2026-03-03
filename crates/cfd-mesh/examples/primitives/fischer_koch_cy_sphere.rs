//! Export: Fischer-Koch C(Y) Sphere → `outputs/primitives/fischer_koch_cy_sphere.stl`
//!
//! Builds a Fischer-Koch C(Y) TPMS lattice mid-surface clipped to a sphere of
//! radius 5 mm, and exports binary STL.
//!
//! The C(Y) (Fischer & Koch 1989):
//! `F = cos(2kx)·sin(ky)·cos(kz) + cos(kx)·cos(2ky)·sin(kz) + sin(kx)·cos(ky)·cos(2kz) = 0`
//!
//! Space group Ia3̄d (#230), genus 6 per unit cell.  Dual Y-shaped triple
//! junctions support triply-connected microfluidic routing and drug delivery.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_fischer_koch_cy_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{FischerKochCySphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Fischer-Koch C(Y) Sphere (TPMS, Marching Cubes)");
    println!("=================================================================");

    let builder = FischerKochCySphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : C(Y) — Fischer & Koch (1989)");
    println!("  Space group: Ia3̄d (#230), genus 6, dual Y-junction labyrinths");
    println!("  Radius     : {} mm", builder.radius);
    println!("  Period     : {} mm  (lattice unit-cell size)", builder.period);
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
    let stl_path = out_dir.join("fischer_koch_cy_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
