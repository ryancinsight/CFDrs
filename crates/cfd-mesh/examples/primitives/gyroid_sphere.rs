//! Export: Gyroid Sphere → `outputs/primitives/gyroid_sphere.stl`
//!
//! Builds a gyroid TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! validates mesh statistics, and exports binary STL.
//!
//! The gyroid is defined by the Schwartz implicit surface:
//!   F(x,y,z) = sin(kx)cos(ky) + sin(ky)cos(kz) + sin(kz)cos(kx) = 0
//!
//! Marching cubes extracts the iso-surface on a voxel grid; only triangles
//! whose centroid is inside the clip sphere are emitted.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_gyroid_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{GyroidSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Gyroid Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = GyroidSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Radius     : {} mm", builder.radius);
    println!(
        "  Period     : {} mm  (lattice unit-cell size)",
        builder.period
    );
    println!(
        "  Resolution : {}  (voxels per axis → {}³ cubes)",
        builder.resolution, builder.resolution
    );
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
    let stl_path = out_dir.join("gyroid_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
