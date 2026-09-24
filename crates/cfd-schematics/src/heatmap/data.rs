use super::svg::build_svg;
use std::path::Path;

/// Data for a single candidate to overlay on the well-plate diagram.
#[derive(Debug, Clone)]
pub struct CandidateZoneData {
    /// Short display label, e.g. `"CCT-lv3"`.
    pub label: String,
    /// Cancer-targeted cavitation score [0–1].
    pub cancer_cav: f64,
    /// Lysis risk index [0–1].
    pub lysis_risk: f64,
    /// Therapy channel / separation efficiency metric [0–1].
    pub therapy_frac: f64,
}

/// Write an SVG 96-well plate diagram to `output_path`.
///
/// Up to 5 `top_candidates` are overlaid as colour bars inside the treatment
/// zone.  Bar colour encodes `cancer_cav` (yellow → red gradient).
///
/// # Errors
///
/// Returns a boxed error if the output file cannot be written.
pub fn write_well_plate_diagram_svg(
    top_candidates: &[CandidateZoneData],
    output_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let svg = build_svg(top_candidates);
    std::fs::write(output_path, svg.as_bytes())?;
    Ok(())
}

