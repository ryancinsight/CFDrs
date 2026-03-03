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

use crate::optimizer::RankedDesign;
use crate::scoring::OptimMode;

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
/// Calls `cfd-schematics`'s `create_geometry` + `plot_geometry` to produce a
/// spatially-laid-out channel plan with x/y mirror symmetry matching the
/// existing `cfd-schematics` example outputs.
///
/// # Errors
/// Returns an error if geometry generation or file writing fails.
pub fn save_schematic_svg(
    candidate: &crate::design::DesignCandidate,
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let system = candidate.to_channel_system();
    let path_str = path.to_str().ok_or("path contains invalid UTF-8")?;
    cfd_schematics::plot_geometry(&system, path_str)?;
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
