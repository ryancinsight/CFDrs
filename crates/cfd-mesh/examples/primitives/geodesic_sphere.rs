//! Export: GeodesicSphere → outputs/primitives/geodesic_sphere.stl
use cfd_mesh::application::watertight::check::check_watertight;
use cfd_mesh::domain::geometry::primitives::PrimitiveMesh;
use cfd_mesh::infrastructure::io::stl;
use cfd_mesh::infrastructure::storage::edge_store::EdgeStore;
use cfd_mesh::GeodesicSphere;
use std::fs;
use std::io::BufWriter;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = GeodesicSphere::default().build()?;
    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
    println!(
        "GeodesicSphere V={} F={} watertight={} vol={:.4} chi={:?}",
        mesh.vertices.len(),
        mesh.faces.len(),
        report.is_watertight,
        report.signed_volume,
        report.euler_characteristic
    );
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("outputs/primitives");
    fs::create_dir_all(&out)?;
    stl::write_binary_stl(
        &mut BufWriter::new(fs::File::create(out.join("geodesic_sphere.stl"))?),
        &mesh.vertices,
        &mesh.faces,
    )?;
    Ok(())
}
