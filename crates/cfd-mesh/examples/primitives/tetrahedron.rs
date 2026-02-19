//! Export: Regular Tetrahedron → `outputs/primitives/tetrahedron.stl`
//!
//! Builds a unit regular tetrahedron using [`cfd_mesh::Tetrahedron`],
//! validates watertightness, and exports binary STL.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_tetrahedron
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::Tetrahedron;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Regular Tetrahedron  (via Tetrahedron builder)");
    println!("=================================================================");

    let mesh = Tetrahedron { radius: 1.0 }.build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Watertight: {}", report.is_watertight);
    println!("  Volume    : {:.6} mm³  (expected ≈ 2.667)", report.signed_volume);
    println!("  Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("tetrahedron.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
