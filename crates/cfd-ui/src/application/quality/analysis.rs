//! Mesh quality analysis — wraps cfd-mesh quality and normal analysis.

use cfd_mesh::application::quality::normals::NormalAnalysis;
use cfd_mesh::application::quality::validation::QualityReport;
use cfd_mesh::IndexedMesh;

/// Consolidated mesh analysis report.
pub struct MeshAnalysisReport {
    /// Quality metrics (aspect ratio, skewness, angle statistics).
    pub quality: QualityReport,
    /// Normal consistency analysis.
    pub normals: NormalAnalysis,
    /// Whether the mesh is watertight.
    pub is_watertight: bool,
    /// Number of vertices.
    pub vertex_count: usize,
    /// Number of faces.
    pub face_count: usize,
    /// Total surface area.
    pub surface_area: f64,
    /// Signed volume (positive = outward normals).
    pub signed_volume: f64,
}

/// Run a full mesh analysis, returning a consolidated report.
pub fn analyze_mesh(mesh: &mut IndexedMesh<f64>) -> MeshAnalysisReport {
    let quality = mesh.quality_report();
    let normals = cfd_mesh::analyze_normals(mesh);
    let is_watertight = mesh.is_watertight();
    MeshAnalysisReport {
        quality,
        normals,
        is_watertight,
        vertex_count: mesh.vertex_count(),
        face_count: mesh.face_count(),
        surface_area: mesh.surface_area(),
        signed_volume: mesh.signed_volume(),
    }
}
