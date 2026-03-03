//! Export: Lidinoid Sphere → `outputs/primitives/lidinoid_sphere.stl`
//!
//! Builds a Lidinoid TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! and exports binary STL.
//!
//! The Lidinoid (Lidin & Larsson 1990):
//! `F = ½(sin(2kx)cos(ky)sin(kz) + sin(2ky)cos(kz)sin(kx) + sin(2kz)cos(kx)sin(ky))`
//! `  − ½(cos(2kx)cos(2ky) + cos(2ky)cos(2kz) + cos(2kz)cos(2kx)) + 0.15 = 0`
//!
//! Space group I4₁32 (#214), genus 3 per unit cell.  Chiral surface — same
//! Bonnet family as the Gyroid (association angle ≈ 38.015°).
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_lidinoid_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{LidinoidSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Lidinoid Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = LidinoidSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : Lidinoid (Lidin & Larsson 1990)");
    println!("  Space group: I4₁32 (#214), genus 3, chiral");
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
    let stl_path = out_dir.join("lidinoid_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
