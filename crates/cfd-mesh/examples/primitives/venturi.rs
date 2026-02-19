//! Export: Venturi Channel → `outputs/primitives/venturi.stl`
//!
//! Builds a millifluidic Venturi tube surface mesh using
//! [`cfd_mesh::geometry::venturi::VenturiMeshBuilder`],
//! validates watertightness, and exports binary STL.
//!
//! The Venturi mesh is aligned along the +Z axis with three named boundary
//! regions: wall (`RegionId(0)`), inlet cap (`RegionId(1)`), outlet cap (`RegionId(2)`).
//!
//! Run with:
//! ```sh
//! cargo run -p cfd-mesh --example primitives_venturi
//! ```

use std::fs;
use std::io::BufWriter;

use cfd_mesh::geometry::venturi::VenturiMeshBuilder;
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=================================================================");
    println!("  Primitive: Venturi Channel  (via VenturiMeshBuilder)");
    println!("=================================================================");
    println!("  d_inlet=2.0 mm  d_throat=0.5 mm");
    println!("  l_inlet=2.0  l_conv=1.0  l_throat=1.0  l_div=1.0  l_outlet=2.0 mm");

    let mesh = VenturiMeshBuilder::new(
        2.0_f64, 0.5, 2.0, 1.0, 1.0, 1.0, 2.0,
    )
    .with_resolution(12, 8)
    .with_circular(true)
    .build_surface()
    .map_err(|e| format!("VenturiMeshBuilder: {e}"))?;

    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);

    println!("  Vertices  : {}", mesh.vertices.len());
    println!("  Faces     : {}", mesh.faces.len());
    println!("  Closed    : {} ({} boundary, {} non-manifold)",
        report.is_closed, report.boundary_edge_count, report.non_manifold_edge_count);
    println!("  Oriented  : {}", report.orientation_consistent);
    println!("  Watertight: {}", report.is_watertight);
    println!("  Volume    : {:.6} mm³", report.signed_volume);
    println!("  Euler χ   : {:?}  (expected 2)", report.euler_characteristic);

    let crate_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let out_dir = crate_dir.join("outputs").join("primitives");
    fs::create_dir_all(&out_dir)?;
    let stl_path = out_dir.join("venturi.stl");

    let file = fs::File::create(&stl_path)?;
    let mut writer = BufWriter::new(file);
    stl::write_binary_stl(&mut writer, &mesh.vertices, &mesh.faces)?;

    println!("  STL       : {}", stl_path.display());
    println!("=================================================================");
    Ok(())
}
