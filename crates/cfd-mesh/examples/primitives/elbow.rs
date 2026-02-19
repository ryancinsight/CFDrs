//! Export: Elbow → outputs/primitives/elbow.stl
use std::fs;
use std::io::BufWriter;
use cfd_mesh::Elbow;
use cfd_mesh::geometry::primitives::PrimitiveMesh;
use cfd_mesh::io::stl;
use cfd_mesh::storage::edge_store::EdgeStore;
use cfd_mesh::watertight::check::check_watertight;
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mesh = Elbow::default().build()?;
    let edges = EdgeStore::from_face_store(&mesh.faces);
    let report = check_watertight(&mesh.vertices, &mesh.faces, &edges);
    println!("Elbow V={} F={} watertight={} vol={:.4} chi={:?}",
        mesh.vertices.len(), mesh.faces.len(), report.is_watertight,
        report.signed_volume, report.euler_characteristic);
    let out = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("outputs/primitives");
    fs::create_dir_all(&out)?;
    stl::write_binary_stl(&mut BufWriter::new(fs::File::create(out.join("elbow.stl"))?),
        &mesh.vertices, &mesh.faces)?;
    Ok(())
}
