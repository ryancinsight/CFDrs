//! Export: Cylinder → `outputs/primitives/cylinder.stl`
//!
//! Builds a closed cylinder (r=1 mm, h=2 mm, 32 segments) using
//! [`cfd_mesh::Cylinder`], validates watertightness, and exports binary STL.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_cylinder
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::Cylinder;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Cylinder r=1 mm, h=2 mm  (via Cylinder builder)");
    println!("=================================================================");

    let mesh = Cylinder {
        radius: 1.0,
        height: 2.0,
        segments: 32,
        ..Cylinder::default()
    }
    .build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
    let expected_vol = PI * 1.0_f64 * 2.0;

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Watertight: {}", report.is_watertight);
    println!(
        "  Volume    : {:.6} mm³  (expected {:.6})",
        report.signed_volume, expected_vol
    );
    println!(
        "  Euler χ   : {:?}  (expected 2)",
        report.euler_characteristic
    );

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("cylinder.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
