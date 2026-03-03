//! STL export command.

use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an `IndexedMesh` to an STL file.
pub fn export_stl(mesh: &IndexedMesh<f64>, path: &Path, binary: bool) -> anyhow::Result<()> {
    let mut file = std::fs::File::create(path)?;
    if binary {
        cfd_mesh::infrastructure::io::stl::write_stl_binary(&mut file, mesh)
            .map_err(|e| anyhow::anyhow!("{e}"))?;
    } else {
        cfd_mesh::infrastructure::io::stl::write_stl_ascii(&mut file, "cfd-ui export", mesh)
            .map_err(|e| anyhow::anyhow!("{e}"))?;
    }
    Ok(())
}
