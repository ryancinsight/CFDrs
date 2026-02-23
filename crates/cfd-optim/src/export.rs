//! JSON and SVG export for top-ranked millifluidic design candidates.
//!
//! # JSON export
//! [`save_top5_json`] serialises a slice of [`RankedDesign`] to a
//! pretty-printed JSON file using `serde_json`.  All metrics and candidate
//! parameters are included.
//!
//! # SVG export
//! [`save_comparison_svg`] produces a horizontal bar chart comparing the
//! scores of the top designs and annotating key metrics (σ, HI/pass,
//! coverage, separation efficiency).  The SVG file can be rendered in any
//! modern browser without external dependencies.

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

// ── SVG export ────────────────────────────────────────────────────────────────

/// Produce an SVG bar-chart comparing the top designs for one optimisation mode.
///
/// The chart includes:
/// - A horizontal score bar (normalised to the highest-scoring candidate).
/// - Annotated key metrics in a table below the bars.
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
    let width: u32 = 900;
    let bar_h: u32 = 48;
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
        "Top designs ranked by composite score — 96-well plate, blood, FDA ≤ 150 Pa",
        &TextStyle::from(("sans-serif", 12).into_font()).color(&RGBColor(80, 80, 80)),
        (20, 46),
    )?;

    // ── Bar area ──────────────────────────────────────────────────────────
    let bar_area = root.margin(header_h as i32, footer_h as i32, 0, 0);
    let max_score = designs
        .iter()
        .map(|d| d.score)
        .fold(0.0_f64, f64::max)
        .max(1e-9);

    // Bar palette — blue gradient by rank
    let bar_colours = [
        RGBColor(41, 128, 185),
        RGBColor(39, 174, 96),
        RGBColor(142, 68, 173),
        RGBColor(243, 156, 18),
        RGBColor(192, 57, 43),
    ];

    let bar_max_px = (width as f64 * 0.48) as i32;

    for (i, d) in designs.iter().enumerate() {
        let y0 = (i as i32) * (bar_h as i32);
        let bar_px = ((d.score / max_score) * bar_max_px as f64) as i32;
        let colour = bar_colours[i % bar_colours.len()];

        // Score bar
        bar_area.draw(&Rectangle::new(
            [(100, y0 + 8), (100 + bar_px, y0 + bar_h as i32 - 8)],
            colour.filled(),
        ))?;

        // Rank label
        bar_area.draw_text(
            &format!("#{}", d.rank),
            &TextStyle::from(("sans-serif", 13).into_font())
                .color(&BLACK),
            (4, y0 + 16),
        )?;

        // Score value
        bar_area.draw_text(
            &format!("{:.4}", d.score),
            &TextStyle::from(("sans-serif", 11).into_font()).color(&BLACK),
            (100 + bar_px + 6, y0 + 16),
        )?;

        // Candidate ID (truncated)
        let id_short = if d.candidate.id.len() > 32 {
            &d.candidate.id[..32]
        } else {
            &d.candidate.id
        };
        bar_area.draw_text(
            id_short,
            &TextStyle::from(("sans-serif", 10).into_font())
                .color(&RGBColor(60, 60, 60)),
            (4, y0 + 32),
        )?;

        // Key metrics annotation (right side)
        let m = &d.metrics;
        let annotation = metric_annotation(d, mode);
        bar_area.draw_text(
            &annotation,
            &TextStyle::from(("sans-serif", 10).into_font())
                .color(&RGBColor(60, 60, 60)),
            (480, y0 + 16),
        )?;

        // FDA badge
        let fda_label = if m.fda_main_compliant { "FDA OK" } else { "FDA !!" };
        let fda_colour = if m.fda_main_compliant {
            RGBColor(39, 174, 96)
        } else {
            RGBColor(192, 57, 43)
        };
        bar_area.draw_text(
            fda_label,
            &TextStyle::from(("sans-serif", 10).into_font()).color(&fda_colour),
            (860, y0 + 16),
        )?;
    }

    // ── Footer ────────────────────────────────────────────────────────────
    root.draw_text(
        &format!(
            "Generated by cfd-optim  |  Plate: ANSI/SLAS 1-2004 96-well  |  {} designs shown",
            designs.len()
        ),
        &TextStyle::from(("sans-serif", 10).into_font())
            .color(&RGBColor(120, 120, 120)),
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
    }
}
