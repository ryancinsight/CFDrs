//! OBJ export command.

use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an `IndexedMesh` to a Wavefront OBJ file.
pub fn export_obj(mesh: &IndexedMesh<f64>, path: &Path) -> anyhow::Result<()> {
    let mut file = std::fs::File::create(path)?;
    cfd_mesh::infrastructure::io::obj::write_obj(&mut file, mesh)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(())
}
