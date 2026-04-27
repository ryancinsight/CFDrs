//! Quality panel — mesh quality metric display.

use crate::application::quality::analysis::MeshAnalysisReport;

/// Format a mesh analysis report as a list of display lines.
#[must_use]
pub fn format_report(report: &MeshAnalysisReport) -> Vec<(String, String)> {
    vec![
        ("Vertices".to_owned(), report.vertex_count.to_string()),
        ("Faces".to_owned(), report.face_count.to_string()),
        (
            "Surface Area".to_owned(),
            format!("{:.6} m^2", report.surface_area),
        ),
        (
            "Signed Volume".to_owned(),
            format!("{:.6} m^3", report.signed_volume),
        ),
        (
            "Watertight".to_owned(),
            if report.is_watertight { "Yes" } else { "No" }.to_owned(),
        ),
    ]
}
