//! Export: UV Sphere → `outputs/primitives/sphere.stl`
//!
//! Builds a UV sphere of radius 1 mm using [`cfd_mesh::UvSphere`],
//! validates watertightness, and exports binary STL.
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_sphere
//! ```

use std::f64::consts::PI;
use std::fs;
use std::io::BufWriter;

use cfd_mesh::UvSphere;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: UV Sphere r=1 mm, 32×16  (via UvSphere builder)");
    println!("=================================================================");

    let sphere = UvSphere {
        radius: 1.0,
        segments: 32,
        stacks: 16,
        ..UvSphere::default()
    };
    let mesh = sphere.build()?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
    let expected_vol = 4.0 * PI / 3.0;

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Watertight: {}", report.is_watertight);
    println!("  Volume    : {:.6} mm³  (expected {:.6})", report.signed_volume, expected_vol);
    println!("  Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("sphere.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
