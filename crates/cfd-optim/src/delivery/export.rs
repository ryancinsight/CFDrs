//! JSON and SVG export for top-ranked millifluidic design candidates.
//!
//! # JSON export
//! [`save_top5_json`] serialises a slice of [`RankedDesign`] to a
//! pretty-printed JSON file using `serde_json`.  All metrics and candidate
//! parameters are included.
//!
//! [`save_all_modes_json`] writes a single combined JSON file containing
//! results from multiple optimisation modes, keyed by mode slug.
//!
//! # SVG export
//! [`save_comparison_svg`] produces a horizontal bar chart comparing the
//! scores of the top designs and annotating key metrics (σ, HI/pass,
//! coverage, separation efficiency).  Bars are plotted on an **absolute
//! [0, 1] scale** so scores remain comparable across modes and runs.  The
//! SVG file can be rendered in any modern browser without external
//! dependencies.

use std::path::Path;

use crate::design::DesignTopology;
use crate::scoring::OptimMode;
use crate::RankedDesign;
use cfd_schematics::visualizations::{plot_blueprint_auto_annotated, RenderConfig};

// ── JSON export ───────────────────────────────────────────────────────────────

/// Serialise a ranked design list to a pretty-printed JSON file.
///
/// # Errors
/// Returns an error if the file cannot be created or serialisation fails.
pub fn save_top5_json(
    designs: &[RankedDesign],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let json = serde_json::to_string_pretty(designs)?;
    std::fs::write(path, json)?;
    Ok(())
}

/// Serialise results from multiple optimisation modes into a single combined
/// JSON file.
///
/// The output is a JSON object of the form:
/// ```json
/// {
///   "cavitation":       [ { "rank": 1, ... }, ... ],
///   "uniform_exposure": [ ... ],
///   ...
/// }
/// ```
///
/// Useful for downstream processing tools and archival, as all modes are
/// available in one round-trip.
///
/// # Errors
/// Returns an error if the file cannot be created or serialisation fails.
pub fn save_all_modes_json(
    mode_results: &[(&str, &[RankedDesign])],
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    use std::collections::BTreeMap;
    // BTreeMap gives stable alphabetical key order in the output JSON.
    let map: BTreeMap<&str, &[RankedDesign]> = mode_results
        .iter()
        .map(|&(slug, designs)| (slug, designs))
        .collect();
    let json = serde_json::to_string_pretty(&map)?;
    std::fs::write(path, json)?;
    Ok(())
}

// ── SVG export ────────────────────────────────────────────────────────────────

/// Produce an SVG bar-chart comparing the top designs for one optimisation mode.
///
/// The chart includes:
/// - A horizontal score bar plotted on an **absolute [0, 1] scale** (not
///   relative to the highest-scoring candidate), allowing scores to be compared
///   across different modes and runs.  A thin gray reference line marks score=1.0.
/// - Annotated key metrics (σ, HI/pass, coverage, separation efficiency).
/// - An **FDA** compliance badge (green OK / red !!) and a **PAI** platelet-
///   activation badge (green OK / amber !!) on the right edge of each row.
///
/// # Errors
/// Returns an error if the file cannot be written.
pub fn save_comparison_svg(
    designs: &[RankedDesign],
    path: &Path,
    mode: OptimMode,
) -> Result<(), Box<dyn std::error::Error>> {
    use plotters::prelude::*;

    if designs.is_empty() {
        return Ok(());
    }

    let mode_name = crate::scoring::score_description(mode);

    // ── Layout constants ──────────────────────────────────────────────────
    let width: u32 = 960;
    let bar_h: u32 = 52;
    let header_h: u32 = 80;
    let footer_h: u32 = 40;
    let n = designs.len() as u32;
    let height: u32 = header_h + n * bar_h + footer_h;

    let root = SVGBackend::new(path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    // ── Title ─────────────────────────────────────────────────────────────
    root.draw_text(
        &format!("cfd-optim  |  {mode_name}"),
        &TextStyle::from(("sans-serif", 18).into_font()).color(&BLACK),
        (20, 20),
    )?;
    root.draw_text(
        "Top designs — absolute score scale [0, 1] — 96-well plate, blood, FDA ≤ 150 Pa",
        &TextStyle::from(("sans-serif", 12).into_font()).color(&RGBColor(80, 80, 80)),
        (20, 46),
    )?;

    // ── Bar area ──────────────────────────────────────────────────────────
    let bar_area = root.margin(header_h as i32, footer_h as i32, 0, 0);

    // Absolute [0, 1] scale: bar_max_px corresponds to score = 1.0
    let bar_max_px: i32 = (f64::from(width) * 0.40) as i32; // 384 px for width=960
    let bar_x0: i32 = 100;

    // Bar palette — distinct colours by rank
    let bar_colours = [
        RGBColor(41, 128, 185), // blue
        RGBColor(39, 174, 96),  // green
        RGBColor(142, 68, 173), // purple
        RGBColor(243, 156, 18), // amber
        RGBColor(192, 57, 43),  // red
    ];

    // Reference line at score=1.0 (thin gray vertical line across all rows)
    let total_bar_h = (n as i32) * (bar_h as i32);
    bar_area.draw(&PathElement::new(
        vec![(bar_x0 + bar_max_px, 0), (bar_x0 + bar_max_px, total_bar_h)],
        RGBColor(200, 200, 200).stroke_width(1),
    ))?;

    for (i, d) in designs.iter().enumerate() {
        let y0 = (i as i32) * (bar_h as i32);
        // Clamp to bar_max_px in case score slightly exceeds 1.0 due to synergy bonuses
        let bar_px = (d.score.min(1.0) * f64::from(bar_max_px)) as i32;
        let colour = bar_colours[i % bar_colours.len()];

        // Score bar
        bar_area.draw(&Rectangle::new(
            [(bar_x0, y0 + 8), (bar_x0 + bar_px, y0 + bar_h as i32 - 10)],
            colour.filled(),
        ))?;

        // Rank label
        bar_area.draw_text(
            &format!("#{}", d.rank),
            &TextStyle::from(("sans-serif", 13).into_font()).color(&BLACK),
            (4, y0 + 16),
        )?;

        // Score value (after bar, or at bar_max+6 for near-full bars)
        let score_x = (bar_x0 + bar_px + 6).min(bar_x0 + bar_max_px + 6);
        bar_area.draw_text(
            &format!("{:.4}", d.score),
            &TextStyle::from(("sans-serif", 11).into_font()).color(&BLACK),
            (score_x, y0 + 16),
        )?;

        // Candidate ID (truncated to 30 chars)
        let id_short = if d.candidate.id.len() > 30 {
            &d.candidate.id[..30]
        } else {
            &d.candidate.id
        };
        bar_area.draw_text(
            id_short,
            &TextStyle::from(("sans-serif", 10).into_font()).color(&RGBColor(60, 60, 60)),
            (4, y0 + 34),
        )?;

        // Key metrics annotation (right-of-bar area)
        let annotation = metric_annotation(d, mode);
        bar_area.draw_text(
            &annotation,
            &TextStyle::from(("sans-serif", 10).into_font()).color(&RGBColor(60, 60, 60)),
            (bar_x0 + bar_max_px + 60, y0 + 16),
        )?;

        // FDA badge
        let m = &d.metrics;
        let fda_label = if m.fda_main_compliant {
            "FDA OK"
        } else {
            "FDA !!"
        };
        let fda_colour = if m.fda_main_compliant {
            RGBColor(39, 174, 96)
        } else {
            RGBColor(192, 57, 43)
        };
        bar_area.draw_text(
            fda_label,
            &TextStyle::from(("sans-serif", 10).into_font()).color(&fda_colour),
            (bar_x0 + bar_max_px + 60, y0 + 34),
        )?;

        // PAI (platelet activation) badge
        let pai_ok = m.platelet_activation_index <= crate::constraints::PAI_PASS_LIMIT;
        let pai_label = if pai_ok { "PAI OK" } else { "PAI !!" };
        let pai_colour = if pai_ok {
            RGBColor(39, 174, 96)
        } else {
            RGBColor(243, 156, 18) // amber warning — PAI overrun is less critical than FDA
        };
        bar_area.draw_text(
            pai_label,
            &TextStyle::from(("sans-serif", 10).into_font()).color(&pai_colour),
            (bar_x0 + bar_max_px + 116, y0 + 34),
        )?;
    }

    // ── Footer ────────────────────────────────────────────────────────────
    root.draw_text(
        &format!(
            "Generated by cfd-optim  |  Plate: ANSI/SLAS 1-2004 96-well  |  {} designs shown  |  bars: absolute score [0,1]",
            designs.len()
        ),
        &TextStyle::from(("sans-serif", 10).into_font()).color(&RGBColor(120, 120, 120)),
        (20, (height - 18) as i32),
    )?;

    root.present()?;
    Ok(())
}

// ── Schematic SVG export ──────────────────────────────────────────────────────

/// Save a 2D channel schematic SVG for a single design candidate.
///
/// The schematic is rendered at the ANSI/SLAS 96-well plate footprint
/// (127.76 × 85.47 mm) with channel geometry positioned within the
/// treatment zone.
///
/// # Errors
/// Returns an error if geometry generation or file writing fails.
pub fn save_schematic_svg(
    candidate: &crate::design::DesignCandidate,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let blueprint = candidate.to_blueprint();
    let path_str = path.to_str().ok_or("path contains invalid UTF-8")?;
    let config = RenderConfig::well_plate_96_report_annotated();
    plot_blueprint_auto_annotated(&blueprint, path_str, &config)?;
    Ok(())
}

/// Build a compact metric annotation string suited to a given optimisation mode.
fn metric_annotation(d: &RankedDesign, mode: OptimMode) -> String {
    let m = &d.metrics;
    match mode {
        OptimMode::SdtCavitation => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "{}  τ_throat={:.0}Pa  HI={:.1e}  cov={:.0}%  ΔP={:.1}kPa",
                sigma,
                m.throat_shear_pa,
                m.hemolysis_index_per_pass,
                m.well_coverage_fraction * 100.0,
                m.total_pressure_drop_pa * 1e-3,
            )
        }
        OptimMode::UniformExposure => {
            format!(
                "unif={:.3}  τ={:.1}Pa  cov={:.0}%  t_res={:.2}s  path={:.1}mm",
                m.flow_uniformity,
                m.max_main_channel_shear_pa,
                m.well_coverage_fraction * 100.0,
                m.mean_residence_time_s,
                m.total_path_length_mm,
            )
        }
        OptimMode::SelectiveAcousticTherapy => {
            format!(
                "3pop={:.3}  ctr={:.0}%  wbc_ctr={:.0}%  rbc_peri={:.0}%  t_res={:.2}s",
                m.three_pop_sep_efficiency,
                m.cancer_center_fraction * 100.0,
                m.wbc_center_fraction * 100.0,
                m.rbc_peripheral_fraction_three_pop * 100.0,
                m.mean_residence_time_s,
            )
        }
        OptimMode::CellSeparation => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "sep={:.3}  {}  HI={:.1e}  cancer_ctr={:.0}%  rbc_peri={:.0}%",
                m.cell_separation_efficiency,
                sigma,
                m.hemolysis_index_per_pass,
                m.cancer_center_fraction * 100.0,
                m.rbc_peripheral_fraction * 100.0,
            )
        }
        OptimMode::ThreePopSeparation => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "3pop={:.3}  {}  wbc_ctr={:.0}%  rbc_peri={:.0}%  HI={:.1e}",
                m.three_pop_sep_efficiency,
                sigma,
                m.wbc_center_fraction * 100.0,
                m.rbc_peripheral_fraction_three_pop * 100.0,
                m.hemolysis_index_per_pass,
            )
        }
        OptimMode::Combined { .. } => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "{}  unif={:.3}  HI={:.1e}  cov={:.0}%",
                sigma,
                m.flow_uniformity,
                m.hemolysis_index_per_pass,
                m.well_coverage_fraction * 100.0,
            )
        }
        OptimMode::SdtTherapy => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "3pop={:.3}  {}  HI={:.1e}  unif={:.3}  cov={:.0}%",
                m.three_pop_sep_efficiency,
                sigma,
                m.hemolysis_index_per_pass,
                m.flow_uniformity,
                m.well_coverage_fraction * 100.0,
            )
        }
        OptimMode::HydrodynamicCavitationSDT => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "dose={:.0}%  {}  rbc@V={:.0}%  3pop={:.3}  HI={:.1e}",
                m.cancer_dose_fraction * 100.0,
                sigma,
                m.rbc_venturi_exposure_fraction * 100.0,
                m.three_pop_sep_efficiency,
                m.hemolysis_index_per_pass,
            )
        }
        OptimMode::PediatricLeukapheresis { patient_weight_kg } => {
            let bv_ml = patient_weight_kg * 85.0;
            format!(
                "wbc_rec={:.0}%  rbc_pass={:.0}%  purity={:.0}%  ECV={:.1}mL/{:.0}mL",
                m.wbc_recovery * 100.0,
                m.rbc_pass_fraction * 100.0,
                m.wbc_purity * 100.0,
                m.total_ecv_ml,
                bv_ml * 0.10,
            )
        }
        OptimMode::CombinedSdtLeukapheresis {
            patient_weight_kg, ..
        } => {
            let bv_ml = patient_weight_kg * 85.0;
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "{}  dose={:.0}%  wbc_rec={:.0}%  ECV={:.1}mL/{:.0}mL  HI={:.1e}  HI15p={:.2}%",
                sigma,
                m.cancer_dose_fraction * 100.0,
                m.wbc_recovery * 100.0,
                m.total_ecv_ml,
                bv_ml * 0.10,
                m.hemolysis_index_per_pass,
                m.projected_hemolysis_15min_pediatric_3kg * 100.0,
            )
        }
        OptimMode::RbcProtectedSdt => {
            let sigma = if m.cavitation_number.is_finite() {
                format!("σ={:.3}", m.cavitation_number)
            } else {
                "σ=∞".to_owned()
            };
            format!(
                "{}  tw={:.3}  lysis={:.2e}  HI={:.1e}  HI15p={:.2}%  fda={}",
                sigma,
                m.therapeutic_window_score,
                m.lysis_risk_index,
                m.hemolysis_index_per_pass,
                m.projected_hemolysis_15min_pediatric_3kg * 100.0,
                if m.fda_overall_compliant {
                    "PASS"
                } else {
                    "FAIL"
                },
            )
        }
    }
}

// ── Annotated primitive selective SVG flow diagram ──────────────────────────

/// Generate an annotated primitive selective topology SVG flow diagram.
///
/// Produces a purpose-built flow diagram showing:
/// - Network topology (pre-trifurcation stages, terminal tri/bi junctions)
/// - Venturi placement on the treatment arm (highlighted)
/// - Per-arm cell population fractions (cancer, WBC, RBC)
/// - Per-channel hemolysis contributions
/// - Flow fractions per arm
///
/// # Errors
/// Returns an error if the file cannot be written.
pub fn save_annotated_selective_svg(
    design: &RankedDesign,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let m = &design.metrics;
    let c = &design.candidate;

    let (topology_label, stage_lines) = match c.topology {
        DesignTopology::PrimitiveSelectiveTree { sequence } => {
            let label = format!("Primitive selective split tree ({})", sequence.label());
            let lines = primitive_stage_lines(sequence);
            (label, lines)
        }
        _ => (format!("{:?}", c.topology), Vec::new()),
    };

    // ANSI/SLAS 96-well plate footprint at 10× scale.
    let width: usize = 1278;
    let row_height = 32;
    let header_height = 120;
    let body_height = stage_lines.len() * row_height + 400;
    let height = header_height + body_height;

    let mut svg = String::with_capacity(8192);
    svg.push_str(&format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">"#
    ));
    svg.push_str(
        r##"<style>
  text { font-family: 'Segoe UI', sans-serif; }
  .title { font-size: 16px; font-weight: bold; fill: #1a1a2e; }
  .subtitle { font-size: 12px; fill: #555; }
  .label { font-size: 11px; fill: #333; }
  .value { font-size: 11px; fill: #0066cc; font-weight: 600; }
  .venturi { font-size: 11px; fill: #cc3300; font-weight: bold; }
  .safe { fill: #228B22; }
  .warn { fill: #cc6600; }
  .stage-box { fill: #f0f4f8; stroke: #aab; rx: 4; ry: 4; }
  .venturi-box { fill: #fff0e0; stroke: #cc6600; rx: 4; ry: 4; stroke-width: 2; }
  .metric-box { fill: #e8f4e8; stroke: #6b9; rx: 4; ry: 4; }
  .arrow { stroke: #666; fill: none; stroke-width: 1.5; marker-end: url(#arrowhead); }
</style>"##,
    );
    svg.push_str(r##"<defs><marker id="arrowhead" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#666"/></marker></defs>"##);

    // Header
    svg.push_str(&format!(
        r#"<text x="20" y="30" class="title">Annotated Primitive Selective Flow Diagram — {}</text>"#,
        design.candidate.id
    ));
    svg.push_str(&format!(
        r#"<text x="20" y="50" class="subtitle">Topology: {topology_label}   Score: {:.4}   Rank: #{}</text>"#,
        design.score, design.rank
    ));

    let sigma_str = if m.cavitation_number.is_finite() {
        format!("σ={:.3}", m.cavitation_number)
    } else {
        "σ=∞".to_owned()
    };
    svg.push_str(&format!(
        r#"<text x="20" y="70" class="subtitle">{sigma_str}   HI/pass={:.2e}   FDA: {}   ΔP={:.1} kPa</text>"#,
        m.hemolysis_index_per_pass,
        if m.fda_overall_compliant { "PASS" } else { "FAIL" },
        m.total_pressure_drop_pa * 1e-3,
    ));

    // Network stage diagram
    let mut y = header_height as f64;
    let x_left = 60.0;

    // Inlet
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="200" height="28" class="stage-box"/>"#,
        x_left,
        y - 14.0
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="label">INLET  Q={:.1} mL/min  HCT={:.1}%</text>"#,
        x_left + 8.0,
        y + 4.0,
        m.flow_rate_ml_min,
        c.feed_hematocrit * 100.0,
    ));
    y += 40.0;

    // Stage lines
    for (i, line) in stage_lines.iter().enumerate() {
        svg.push_str(&format!(
            r#"<line x1="{}" y1="{}" x2="{}" y2="{}" class="arrow"/>"#,
            x_left + 100.0,
            y - 30.0,
            x_left + 100.0,
            y - 14.0,
        ));
        svg.push_str(&format!(
            r#"<rect x="{}" y="{}" width="500" height="28" class="stage-box"/>"#,
            x_left,
            y - 14.0
        ));
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" class="label">Stage {}: {line}</text>"#,
            x_left + 8.0,
            y + 4.0,
            i + 1,
        ));
        y += row_height as f64 + 8.0;
    }

    // Venturi throat
    svg.push_str(&format!(
        r#"<line x1="{}" y1="{}" x2="{}" y2="{}" class="arrow"/>"#,
        x_left + 100.0,
        y - 30.0,
        x_left + 100.0,
        y - 14.0,
    ));
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="500" height="28" class="venturi-box"/>"#,
        x_left,
        y - 14.0
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="venturi">▶ VENTURI THROAT  d={:.0}µm  L={:.0}µm  τ={:.0} Pa  t={:.2e} s</text>"#,
        x_left + 8.0,
        y + 4.0,
        c.throat_diameter_m * 1e6,
        c.throat_length_m * 1e6,
        m.throat_shear_pa,
        m.throat_transit_time_s,
    ));
    y += 50.0;

    // Cell population summary
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="820" height="120" class="metric-box"/>"#,
        x_left - 20.0,
        y - 10.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="title">Cell Population Fractions at Venturi</text>"#,
        x_left,
        y + 14.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Cancer center fraction: {:.1}%</text>"#,
        x_left,
        y + 34.0,
        m.cancer_center_fraction * 100.0,
    ));
    svg.push_str(&format!(
        r#"<text x="350" y="{}" class="value">RBC peripheral fraction: {:.1}%</text>"#,
        y + 34.0,
        m.rbc_peripheral_fraction * 100.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">WBC center fraction: {:.1}%</text>"#,
        x_left,
        y + 54.0,
        m.wbc_center_fraction * 100.0,
    ));
    svg.push_str(&format!(
        r#"<text x="350" y="{}" class="value">Local HCT at venturi: {:.1}%</text>"#,
        y + 54.0,
        m.local_hematocrit_venturi * 100.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Venturi flow fraction: {:.1}%</text>"#,
        x_left,
        y + 74.0,
        m.venturi_flow_fraction * 100.0,
    ));
    svg.push_str(&format!(
        r#"<text x="350" y="{}" class="value">RBC venturi exposure: {:.1}%</text>"#,
        y + 74.0,
        m.rbc_venturi_exposure_fraction * 100.0,
    ));
    y += 140.0;

    // Hemolysis decomposition
    svg.push_str(&format!(
        r#"<rect x="{}" y="{}" width="820" height="100" class="metric-box"/>"#,
        x_left - 20.0,
        y - 10.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="title">Per-Channel Hemolysis Decomposition</text>"#,
        x_left,
        y + 14.0,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Treatment channel HI: {:.2e}</text>"#,
        x_left,
        y + 34.0,
        m.treatment_channel_hi,
    ));
    svg.push_str(&format!(
        r#"<text x="450" y="{}" class="value">Bulk HI (uncorrected): {:.2e}</text>"#,
        y + 34.0,
        m.bulk_hemolysis_index_per_pass,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Bypass mean HI: {:.2e}</text>"#,
        x_left,
        y + 54.0,
        m.bypass_channel_hi_mean,
    ));
    svg.push_str(&format!(
        r#"<text x="450" y="{}" class="value">Final HI (corrected): {:.2e}</text>"#,
        y + 54.0,
        m.hemolysis_index_per_pass,
    ));
    svg.push_str(&format!(
        r#"<text x="{}" y="{}" class="value">Bypass max HI: {:.2e}</text>"#,
        x_left,
        y + 74.0,
        m.bypass_channel_hi_max,
    ));
    let hi_reduction_pct = if m.bulk_hemolysis_index_per_pass > 1e-18 {
        (1.0 - m.hemolysis_index_per_pass / m.bulk_hemolysis_index_per_pass) * 100.0
    } else {
        0.0
    };
    let class = if hi_reduction_pct > 0.0 {
        "safe"
    } else {
        "warn"
    };
    svg.push_str(&format!(
        r#"<text x="450" y="{}" class="value {class}">HCT correction: {:.1}% reduction</text>"#,
        y + 74.0,
        hi_reduction_pct,
    ));

    // Primitive selective specific flow fractions
    y += 120.0;
    if matches!(c.topology, DesignTopology::PrimitiveSelectiveTree { .. }) {
        svg.push_str(&format!(
            r#"<rect x="{}" y="{}" width="820" height="60" class="stage-box"/>"#,
            x_left - 20.0,
            y - 10.0,
        ));
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" class="label">Selective model Q_venturi={:.1}%   solved Q_venturi={:.1}%   pretri_qfrac_mean={:.3}   tri_qfrac={:.3}   bi_qfrac={:.3}</text>"#,
            x_left,
            y + 14.0,
            m.cif_model_venturi_flow_fraction * 100.0,
            m.cif_solved_venturi_flow_fraction * 100.0,
            m.cif_pretri_qfrac_mean,
            m.cif_terminal_tri_qfrac,
            m.cif_terminal_bi_qfrac,
        ));
        svg.push_str(&format!(
            r#"<text x="{}" y="{}" class="label">Outlet tail: {:.1} mm   Remerge proximity: {:.3}   SCDI: {:.4}</text>"#,
            x_left,
            y + 36.0,
            m.cif_outlet_tail_length_mm,
            m.cif_remerge_proximity_score,
            m.selective_cavitation_delivery_index,
        ));
    }

    svg.push_str("</svg>");
    std::fs::write(path, svg)?;
    Ok(())
}

fn primitive_stage_lines(sequence: crate::design::PrimitiveSplitSequence) -> Vec<String> {
    let bits = sequence.split_types();
    (0..usize::from(sequence.levels()))
        .map(|idx| {
            if ((bits >> idx) & 1) == 0 {
                format!(
                    "Bifurcation stage {}: treatment branch remains selected, bypass branch merges later",
                    idx + 1
                )
            } else {
                format!(
                    "Trifurcation stage {}: center treatment band continues, peripheral bands bypass therapy",
                    idx + 1
                )
            }
        })
        .collect()
}
