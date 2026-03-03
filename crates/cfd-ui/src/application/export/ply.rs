//! PLY export command.

use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an `IndexedMesh` to a Stanford PLY file.
pub fn export_ply(mesh: &IndexedMesh<f64>, path: &Path) -> anyhow::Result<()> {
    let mut file = std::fs::File::create(path)?;
    cfd_mesh::infrastructure::io::ply::write_ply(&mut file, mesh)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(())
}
