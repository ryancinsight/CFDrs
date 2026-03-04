//! Design Option 1 (ESDT) CIF schematic SVG generator.
//!
//! Produces an annotated CIF topology diagram for the "External Sonodynamic
//! Therapy" design option: [`IncrementalFiltrationTriBiSeparator`] with
//! **no venturi** — ultrasound-only treatment.
//!
//! Shows channel widths, flow paths, Zweifach-Fung routing, and the
//! ultrasound treatment zone on the center arm.

use std::path::Path;

/// Parameters for the Design Option 1 CIF schematic.
pub struct DesignOption1Params {
    /// Number of pre-trifurcation stages (typically 2).
    pub n_pretri: u8,
    /// Parent channel width \[mm\].
    pub parent_width_mm: f64,
    /// Center-arm fraction at pre-trifurcation stages.
    pub pretri_center_frac: f64,
    /// Center-arm fraction at the terminal trifurcation.
    pub terminal_tri_center_frac: f64,
    /// Treatment-arm fraction at the terminal bifurcation.
    pub terminal_bi_treat_frac: f64,
    /// Feed hematocrit.
    pub feed_hematocrit: f64,
    /// Total flow rate \[mL/min\].
    pub flow_rate_ml_min: f64,
}

impl Default for DesignOption1Params {
    fn default() -> Self {
        Self {
            n_pretri: 2,
            parent_width_mm: 6.0,
            pretri_center_frac: 0.50,
            terminal_tri_center_frac: 0.53,
            terminal_bi_treat_frac: 0.80,
            feed_hematocrit: 0.45,
            flow_rate_ml_min: 100.0,
        }
    }
}

/// Generate a Design Option 1 (ESDT, no-venturi) CIF schematic SVG.
///
/// The diagram illustrates the full CIF staging architecture:
/// 1. `n_pretri` pre-trifurcation stages with asymmetric center bias
/// 2. Terminal trifurcation for final cell enrichment
/// 3. Terminal bifurcation routing the treatment arm to the
///    **ultrasound-only treatment zone** (no hydrodynamic cavitation)
/// 4. Bypass arms merging to external outlet
///
/// Each stage is annotated with channel widths and flow fractions.
///
/// # Errors
/// Returns an error if the file cannot be written.
pub fn save_design_option1_svg(
    params: &DesignOption1Params,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    // ANSI/SLAS 96-well plate footprint: 127.76 × 85.47 mm.
    // Pixel dimensions at 10× scale preserve the plate aspect ratio.
    const PLATE_PX_W: usize = 1278;
    const PLATE_PX_H: usize = 855;

    let width = PLATE_PX_W;
    let height = PLATE_PX_H;

    let mut svg = String::with_capacity(8192);
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#
    ));
    svg.push_str(STYLE);

    // Header
    svg.push_str(&format!(
        r#"<text x="20" y="30" class="title">Design Option 1 — ESDT CIF Schematic (Ultrasound-Only SDT)</text>"#
    ));
    svg.push_str(&format!(
        r#"<text x="20" y="52" class="subtitle">IncrementalFiltrationTriBiSeparator  n_pretri={}  Q={:.0} mL/min  HCT={:.0}%  No Venturi</text>"#,
        params.n_pretri,
        params.flow_rate_ml_min,
        params.feed_hematocrit * 100.0,
    ));

    let mut y: f64 = 80.0;
    let x_left = 60.0;
    let box_w = 780.0;
    let arrow_x = 200.0;

    // Cumulative flow fraction tracking: center carries `q`, bypass collects rest.
    let mut q_center = 1.0_f64;
    let mut current_width_mm = params.parent_width_mm;

    // Inlet box
    draw_stage(
        &mut svg,
        x_left,
        y,
        box_w,
        &format!(
            "INLET — W={:.1} mm — Q={:.0} mL/min — HCT={:.0}%",
            params.parent_width_mm,
            params.flow_rate_ml_min,
            params.feed_hematocrit * 100.0,
        ),
        "stage-box",
    );
    y += 50.0;

    // Pre-trifurcation stages
    for i in 0..params.n_pretri {
        draw_arrow(&mut svg, arrow_x, y - 20.0, y);
        draw_junction(&mut svg, arrow_x, y, i as usize, 1.3);

        let center_w = current_width_mm * params.pretri_center_frac;
        let periph_w = current_width_mm * (1.0 - params.pretri_center_frac) / 2.0;
        let q_periph = q_center * (1.0 - params.pretri_center_frac);
        q_center *= params.pretri_center_frac;

        let label = format!(
            "Pre-tri stage {} — Zweifach-Fung split: center W={:.2} mm ({:.0}%), \
             periphery 2×{:.2} mm → bypass  |  Q_center={:.1}%",
            i + 1,
            center_w,
            params.pretri_center_frac * 100.0,
            periph_w,
            q_center * 100.0,
        );
        draw_stage(&mut svg, x_left, y, box_w, &label, "stage-box");

        // Bypass annotation (right side)
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" class="bypass">→ bypass ({:.1}% flow, RBC-enriched)</text>"#,
            x_left + box_w + 10.0,
            y + 18.0,
            q_periph * 100.0,
        ));

        current_width_mm = center_w;
        y += 80.0;
    }

    // Terminal trifurcation
    draw_arrow(&mut svg, arrow_x, y - 20.0, y);
    draw_junction(&mut svg, arrow_x, y, params.n_pretri as usize, 1.3);
    let tri_center_w = current_width_mm * params.terminal_tri_center_frac;
    let tri_periph_w = current_width_mm * (1.0 - params.terminal_tri_center_frac) / 2.0;
    let q_tri_periph = q_center * (1.0 - params.terminal_tri_center_frac);
    q_center *= params.terminal_tri_center_frac;

    draw_stage(
        &mut svg,
        x_left,
        y,
        box_w,
        &format!(
            "Terminal trifurcation — center W={:.2} mm ({:.0}%), \
             periphery 2×{:.2} mm → bypass  |  Q_center={:.1}%",
            tri_center_w,
            params.terminal_tri_center_frac * 100.0,
            tri_periph_w,
            q_center * 100.0,
        ),
        "stage-box",
    );
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="bypass">→ bypass ({:.1}% flow)</text>"#,
        x_left + box_w + 10.0,
        y + 18.0,
        q_tri_periph * 100.0,
    ));
    y += 80.0;

    // Terminal bifurcation
    draw_arrow(&mut svg, arrow_x, y - 20.0, y);
    draw_junction(&mut svg, arrow_x, y, params.n_pretri as usize + 1, 1.5);
    let treat_w = tri_center_w * params.terminal_bi_treat_frac;
    let waste_w = tri_center_w * (1.0 - params.terminal_bi_treat_frac);
    let q_waste = q_center * (1.0 - params.terminal_bi_treat_frac);
    q_center *= params.terminal_bi_treat_frac;

    draw_stage(
        &mut svg,
        x_left,
        y,
        box_w,
        &format!(
            "Terminal bifurcation — treatment arm W={:.2} mm ({:.0}%), \
             waste arm W={:.2} mm → bypass  |  Q_treat={:.1}%",
            treat_w,
            params.terminal_bi_treat_frac * 100.0,
            waste_w,
            q_center * 100.0,
        ),
        "stage-box",
    );
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="bypass">→ bypass ({:.1}% flow)</text>"#,
        x_left + box_w + 10.0,
        y + 18.0,
        q_waste * 100.0,
    ));
    y += 80.0;

    // Ultrasound treatment zone (no venturi)
    draw_arrow(&mut svg, arrow_x, y - 20.0, y);
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="{}" height="50" class="us-zone"/>"#,
        x_left,
        y - 4.0,
        box_w,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="us-label">◈ ULTRASOUND TREATMENT ZONE — \
        No Venturi — Cancer-enriched center stream — Q={:.1}%  W={:.2} mm</text>"#,
        x_left + 12.0,
        y + 18.0,
        q_center * 100.0,
        treat_w,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="subtitle">Selective cancer-cell exposure via external ultrasound transducer</text>"#,
        x_left + 12.0,
        y + 38.0,
    ));
    y += 80.0;

    // Outlet merge
    draw_arrow(&mut svg, arrow_x, y - 20.0, y);
    draw_stage(
        &mut svg,
        x_left,
        y,
        box_w,
        "OUTLET — All bypass arms merge → single external port",
        "stage-box",
    );
    y += 60.0;

    // Summary box
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="{}" height="80" class="summary-box"/>"#,
        x_left - 10.0,
        y,
        box_w + 20.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="title">Design Option 1 Summary</text>"#,
        x_left + 4.0,
        y + 20.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Treatment flow: {:.1}% of inlet  |  Treatment channel: {:.2} mm  |  {} CIF stages total</text>"#,
        x_left + 4.0,
        y + 42.0,
        q_center * 100.0,
        treat_w,
        params.n_pretri as usize + 2, // pretri + terminal tri + terminal bi
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Mechanism: Ultrasound-only SDT (no hydrodynamic cavitation)  |  Zweifach-Fung routing: cancer → center, RBC → periphery</text>"#,
        x_left + 4.0,
        y + 62.0,
    ));

    svg.push_str("</svg>");
    std::fs::write(path, svg)?;
    Ok(())
}

fn draw_stage(svg: &mut String, x: f64, y: f64, w: f64, label: &str, class: &str) {
    svg.push_str(&format!(
        r#"<rect x="{x}" y="{}" width="{w}" height="36" class="{class}"/>"#,
        y - 4.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="label">{label}</text>"#,
        x + 10.0,
        y + 18.0,
    ));
}

fn draw_arrow(svg: &mut String, x: f64, y1: f64, y2: f64) {
    svg.push_str(&format!(
        r#"<line x1="{x}" y1="{y1}" x2="{x}" y2="{y2}" class="arrow"/>"#
    ));
}

/// Draw a junction indicator circle + K-factor annotation at a channel intersection.
///
/// # Arguments
/// * `stage` — stage index (0-based), used for unique element IDs
/// * `k_factor` — Idelchik K-factor for this junction type
fn draw_junction(svg: &mut String, x: f64, y: f64, stage: usize, k_factor: f64) {
    // Yellow circle with red stroke marks the physical intersection node.
    svg.push_str(&format!(
        r#"<circle id="jn_{stage}" cx="{x}" cy="{}" r="8" class="jn-marker"/>"#,
        y - 10.0,
    ));
    // K-factor annotation to the right of the junction circle.
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="jn-label">K={k_factor:.1} (Idelchik)</text>"#,
        x + 14.0,
        y - 5.0,
    ));
}

const STYLE: &str = r##"<style>
  text { font-family: 'Segoe UI', sans-serif; }
  .title { font-size: 15px; font-weight: bold; fill: #1a1a2e; }
  .subtitle { font-size: 11px; fill: #555; }
  .label { font-size: 11px; fill: #333; }
  .value { font-size: 11px; fill: #0066cc; font-weight: 600; }
  .bypass { font-size: 10px; fill: #888; font-style: italic; }
  .us-label { font-size: 12px; fill: #7b2d8e; font-weight: bold; }
  .stage-box { fill: #f0f4f8; stroke: #aab; rx: 4; ry: 4; }
  .us-zone { fill: #f3e5f5; stroke: #9c27b0; rx: 6; ry: 6; stroke-width: 2; stroke-dasharray: 6,3; }
  .summary-box { fill: #e8f4e8; stroke: #4caf50; rx: 6; ry: 6; }
  .arrow { stroke: #666; fill: none; stroke-width: 1.5; marker-end: url(#arrowhead); }
  .jn-marker { fill: #fff176; stroke: #e53935; stroke-width: 2; }
  .jn-label { font-size: 9px; fill: #e53935; font-weight: 600; }
</style>
<defs><marker id="arrowhead" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#666"/></marker></defs>"##;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_valid_svg() {
        let dir = std::env::temp_dir().join("cfd_cif_schematic_test");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("design_option1.svg");

        save_design_option1_svg(&DesignOption1Params::default(), &path).unwrap();

        let content = std::fs::read_to_string(&path).unwrap();
        assert!(content.starts_with("<svg"));
        assert!(content.ends_with("</svg>"));
        assert!(content.contains("ULTRASOUND TREATMENT ZONE"));
        assert!(content.contains("No Venturi"));
        assert!(content.contains("Pre-tri stage 1"));
        assert!(content.contains("Pre-tri stage 2"));
        assert!(content.contains("Terminal trifurcation"));
        assert!(content.contains("Terminal bifurcation"));

        std::fs::remove_dir_all(&dir).ok();
    }
}
