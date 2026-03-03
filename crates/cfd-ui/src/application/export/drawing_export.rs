//! Engineering drawing export command.

use crate::domain::drawing::sheet::DrawingSheet;
use cfd_mesh::IndexedMesh;
use std::path::Path;

/// Export an engineering drawing to SVG.
pub fn export_drawing_svg(
    sheet: &DrawingSheet,
    meshes: &[&IndexedMesh<f64>],
    path: &Path,
) -> anyhow::Result<()> {
    let svg_content = crate::infrastructure::drawing::svg_renderer::render_svg(sheet, meshes);
    std::fs::write(path, svg_content)?;
    Ok(())
}
