//! glTF 2.0 Binary (GLB) export command.

use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an `IndexedMesh` to a glTF Binary (`.glb`) file.
pub fn export_glb(mesh: &IndexedMesh<f64>, path: &Path) -> anyhow::Result<()> {
    let mut file = std::fs::File::create(path)?;
    cfd_mesh::infrastructure::io::gltf_export::write_glb(&mut file, mesh)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(())
}
