//! Export: Axis-Aligned Cube → `outputs/primitives/cube.stl`
//!
//! Builds a 2 mm³ cube using [`cfd_mesh::Cube`], validates watertightness,
//! and exports binary STL.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_cube
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::core::scalar::Point3r;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::Cube;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Cube 2×2×2 mm  (via Cube builder)");
    println!("=================================================================");

    let mesh = Cube {
        origin: Point3r::origin(),
        width: 2.0,
        height: 2.0,
        depth: 2.0,
    }
    .build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Watertight: {}", report.is_watertight);
    println!(
        "  Volume    : {:.6} mm³  (expected 8.0)",
        report.signed_volume
    );
    println!(
        "  Euler χ   : {:?}  (expected 2)",
        report.euler_characteristic
    );

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("cube.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
