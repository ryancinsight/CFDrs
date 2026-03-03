//! DXF export command.

use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an `IndexedMesh` to a DXF file with 3DFACE entities.
pub fn export_dxf(mesh: &IndexedMesh<f64>, path: &Path) -> anyhow::Result<()> {
    let mut file = std::fs::File::create(path)?;
    cfd_mesh::infrastructure::io::dxf::write_dxf(&mut file, mesh)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(())
}
