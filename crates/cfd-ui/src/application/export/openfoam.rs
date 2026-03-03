//! OpenFOAM polyMesh export command.

use cfd_mesh::IndexedMesh;
use cfd_mesh::domain::core::index::RegionId;
use cfd_mesh::domain::topology::halfedge::PatchType;
use std::path::Path;

/// Export an `IndexedMesh` to OpenFOAM polyMesh format in the given directory.
///
/// `patches` maps region IDs to patch names and types. Faces whose region ID
/// is not listed are collected into a synthetic `"defaultFaces"` wall patch.
pub fn export_openfoam(
    mesh: &IndexedMesh<f64>,
    dir: &Path,
    patches: &[(RegionId, &str, PatchType)],
) -> anyhow::Result<()> {
    std::fs::create_dir_all(dir)?;
    cfd_mesh::infrastructure::io::openfoam::write_openfoam_polymesh(mesh, dir, patches)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    Ok(())
}
