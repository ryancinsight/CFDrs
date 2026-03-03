//! Export: FRD Sphere → `outputs/primitives/frd_sphere.stl`
//!
//! Builds a Schoen FRD TPMS lattice mid-surface clipped to a sphere of radius 5 mm,
//! and exports binary STL.
//!
//! The FRD (Schoen 1970 — Face-centred Rhombic Dodecahedron):
//! `F = 4·cos(kx)·cos(ky)·cos(kz) − (cos(2kx)·cos(2ky)+cos(2ky)·cos(2kz)+cos(2kz)·cos(2kx)) = 0`
//!
//! FCC topology with two congruent interlocking labyrinths.  Wider channel
//! junctions than Schwarz D: high-flow filtration and fuel cell scaffold use.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_frd_sphere --release
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::domain::geometry::primitives::{FrdSphere, PrimitiveMesh};
use cfd_mesh::infrastructure::io::stl;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: FRD Sphere (TPMS mid-surface, Marching Cubes)");
    println!("=================================================================");

    let builder = FrdSphere {
        radius: 5.0,
        period: 2.5,
        resolution: 80,
        iso_value: 0.0,
    };

    println!("  Surface    : FRD — Schoen F-RD (1970), Face-centred Rhombic Dodecahedron");
    println!("  Topology   : FCC, two congruent interlocking labyrinths");
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
    let stl_path = out_dir.join("frd_sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL        : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
