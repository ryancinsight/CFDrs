//! Export: Split P Sphere → `outputs/primitives/split_p_sphere.stl`
//!
//! Builds a Split P TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! and exports binary STL.
//!
//! The Split P (Fogden & Hyde 1992):
//! `F = 1.1(sin(2kx)sin(kz)cos(ky)+sin(2ky)sin(kx)cos(kz)+sin(2kz)sin(ky)cos(kx))`
//! `  − 0.2(cos(2kx)cos(2ky)+cos(2ky)cos(2kz)+cos(2kz)cos(2kx))`
//! `  − 0.4(cos(2kx)+cos(2ky)+cos(2kz)) = 0`
//!
//! Tetragonal symmetry; rPD Bonnet family.  Anisotropic channel diameters
//! enable directional permeability control in microfluidic membranes.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_split_p_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{SplitPSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Split P Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = SplitPSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : Split P (Fogden & Hyde 1992, rPD Bonnet family)");
    println!("  Symmetry   : Tetragonal, anisotropic channels");
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
    let stl_path = out_dir.join("split_p_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
