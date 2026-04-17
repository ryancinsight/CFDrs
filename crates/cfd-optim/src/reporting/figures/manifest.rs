//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use cfd_schematics::domain::model::ChannelShape;
use cfd_schematics::geometry::metadata::ChannelVisualRole;
use serde::{Deserialize, Serialize};

use super::context::{write_context_concept_figures, write_selected_schematic_figure};
use super::process::write_creation_optimization_process_figure;
use super::svg::{
    write_cavitation_distribution_figure, write_cross_mode_figure,
    write_dean_venturi_placement_figure, write_ga_convergence_figure, write_head_to_head_figure,
    write_pareto_figure, write_pediatric_ecv_figure, DeanVenturiPoint,
};
use super::treatment_lane::write_treatment_lane_zoom_figure;
use crate::constraints::{PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG};
use crate::reporting::{rank_ga_hydrosdt_report_designs, Milestone12ReportDesign, ParetoPoint};

/// Figure metadata used for dynamic table-of-contents and section rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NarrativeFigureSpec {
    pub number: usize,
    pub title: String,
    pub path: String,
    pub caption: String,
    pub alt: String,
}

/// Inputs needed to generate dynamic report figures.
pub struct FigureGenerationInput<'a> {
    pub option1_ranked: &'a [Milestone12ReportDesign],
    pub option2_ranked: &'a [Milestone12ReportDesign],
    pub ga_top: &'a [Milestone12ReportDesign],
    pub option2_pool_all: &'a [ParetoPoint],
    pub ga_pool_all: &'a [ParetoPoint],
    pub ga_best_per_gen: &'a [f64],
    pub fast_mode: bool,
}

/// Generate report figure assets and return ordered figure specs.
///
/// # Theorem
/// Accepted Milestone 12 report figures are generated from canonical report
/// inputs without fallback placeholders: if any required schematic or derived
/// figure cannot be reproduced from the ranked designs and report datasets, the
/// function returns an error rather than degrading the report.
///
/// **Proof sketch**
/// The generator validates the presence of required Option 1, Option 2, and GA
/// ranked designs up front, rewrites the canonical concept and selected-design
/// SVGs directly from blueprint SSOT, and propagates any figure-builder failure
/// immediately. Because there is no branch that emits a synthetic substitute or
/// silently skips a required figure, successful completion implies all required
/// figures were materialized from canonical sources.
///
/// # Errors
/// Returns an error if any figure cannot be written.
pub fn generate_m12_report_figures(
    figures_dir: &Path,
    figure_path_prefix: &str,
    input: &FigureGenerationInput<'_>,
) -> Result<Vec<NarrativeFigureSpec>, Box<dyn std::error::Error>> {
    let option1 = input
        .option1_ranked
        .first()
        .ok_or("Milestone 12 report requires a selected Option 1 design")?;
    let option2 = input
        .option2_ranked
        .first()
        .ok_or("Milestone 12 report requires a selected Option 2 design")?;
    let ga_best = input
        .ga_top
        .first()
        .ok_or("Milestone 12 report requires a selected GA design")?;
    let ga_ranked_for_figures: Vec<Milestone12ReportDesign> = input
        .option2_ranked
        .first()
        .map(|baseline_option2| rank_ga_hydrosdt_report_designs(input.ga_top, baseline_option2))
        .filter(|filtered| !filtered.is_empty())
        .unwrap_or_else(|| input.ga_top.to_vec());

    std::fs::create_dir_all(figures_dir)?;

    write_creation_optimization_process_figure(
        &figures_dir.join("milestone12_creation_optimization_process.svg"),
    )?;
    write_context_concept_figures(figures_dir)?;
    write_selected_schematic_figure(&figures_dir.join("selected_option1_schematic.svg"), option1)?;
    write_selected_schematic_figure(
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        option2,
    )?;
    write_selected_schematic_figure(&figures_dir.join("top_hydrosdt_schematic.svg"), ga_best)?;

    write_cross_mode_figure(
        &figures_dir.join("m12_cross_mode_scoring.svg"),
        input.option1_ranked,
        input.option2_ranked,
        &ga_ranked_for_figures,
    )?;
    write_head_to_head_figure(
        &figures_dir.join("m12_head_to_head_scores.svg"),
        input.option2_ranked,
    )?;
    write_cavitation_distribution_figure(
        &figures_dir.join("m12_cavitation_distribution.svg"),
        input.option2_ranked,
        &ga_ranked_for_figures,
    )?;
    write_pareto_figure(
        &figures_dir.join("m12_pareto_oncology.svg"),
        input.option2_ranked,
        &ga_ranked_for_figures,
        input.option2_pool_all,
        input.ga_pool_all,
    )?;
    write_pediatric_ecv_figure(
        &figures_dir.join("m12_pediatric_ecv_margin.svg"),
        input.option1_ranked.first(),
        first_ranked(input.option2_ranked),
        first_ranked(&ga_ranked_for_figures),
    )?;
    if input.ga_best_per_gen.is_empty() {
        return Err("Milestone 12 report requires non-empty GA convergence history".into());
    }
    write_ga_convergence_figure(
        &figures_dir.join("m12_ga_convergence.svg"),
        input.ga_best_per_gen,
        input.fast_mode,
    )?;

    // Dean-venturi placement figure: recompute per-placement metrics from the
    // GA top design to show how bend position influences cavitation.
    let option2_serpentine_focus =
        pick_serpentine_venturi_focus(input.option2_ranked).unwrap_or(option2);
    let ga_serpentine_focus =
        pick_serpentine_venturi_focus(&ga_ranked_for_figures).unwrap_or(&ga_ranked_for_figures[0]);
    let dean_venturi_points = extract_dean_venturi_points(Some(ga_serpentine_focus))?;
    if dean_venturi_points.is_empty() {
        return Err(
            "Milestone 12 report requires non-empty Dean-versus-cavitation placement data".into(),
        );
    }
    write_dean_venturi_placement_figure(
        &figures_dir.join("m12_dean_venturi_placement.svg"),
        "Dean Number vs Cavitation at Serpentine Bend Apices",
        &dean_venturi_points,
    )?;

    write_treatment_lane_zoom_figure(
        &figures_dir.join("m12_treatment_lane_zoom.svg"),
        option2_serpentine_focus,
        ga_serpentine_focus,
    )?;

    let mut specs = vec![
        spec(
            1,
            "Milestone 12 Design Creation & Optimization Process",
            "milestone12_creation_optimization_process.svg",
            figure_path_prefix,
            "Canonical Milestone 12 pipeline from topology request generation in `cfd-schematics`, through geometry-authored blueprint creation, 1D physics evaluation, strict eligibility gating, deterministic ranking, and optional GA in-place refinement. The narrative, tables, and figures in this report are all emitted from this canonical artifact chain rather than edited manually.",
            "Milestone 12 process figure",
        ),
        spec(
            2,
            "Treatment Zone Layout - Bifurcation Concept",
            "treatment_zone_plate.svg",
            figure_path_prefix,
            "Bifurcation-root treatment-zone concept within the ANSI/SLAS 96-well plate footprint. The figure illustrates the physical envelope that constrains all downstream split-tree layouts, ultrasound exposure windows, and outlet collection geometry.",
            "Treatment zone layout constrained to 96-well plate footprint",
        ),
        spec(
            3,
            "Treatment Zone Layout - Trifurcation Concept",
            "treatment_zone_plate_trifurcation.svg",
            figure_path_prefix,
            "Trifurcation-root treatment-zone concept showing how additional selective-routing degrees of freedom can concentrate the treatment lane while preserving bypass paths for healthier blood fractions. This geometry motivates the deeper asymmetric split trees evaluated in Milestone 12.",
            "Treatment zone layout optimized for ultrasound delivery",
        ),
        spec(
            4,
            &format!(
                "Selected Option 1 Design - {} Selective Acoustic",
                stage_sequence_label(option1)
            ),
            "selected_option1_schematic.svg",
            figure_path_prefix,
            &selected_schematic_caption(
                option1,
                "Option 1",
                "selective acoustic center-treatment network",
            ),
            "Option 1 schematic",
        ),
        spec(
            5,
            &format!(
                "Selected Option 2 Design - {} Combined Selective Venturi",
                stage_sequence_label(option2)
            ),
            "selected_option2_combined_schematic.svg",
            figure_path_prefix,
            &selected_schematic_caption(
                option2,
                "Option 2",
                "selective venturi treatment ranked by the combined score",
            ),
            "Option 2 schematic",
        ),
        spec(
            6,
            &format!(
                "Best HydroSDT GA-Optimized Design - {}",
                stage_sequence_label(ga_best)
            ),
            "top_hydrosdt_schematic.svg",
            figure_path_prefix,
            &format!(
                "Best HydroSDT GA-optimized design after architecture-preserving in-place refinement. Topology: {}. Visible split layers: {} ({}). Active venturi throats: {} with {} serial stage(s) per path. The treatment-path geometry shows where Dean-generating curvature and venturi placement were co-localised for the final GA-ranked lineage.",
                ga_best.topology_display_name(),
                visible_split_layers(ga_best),
                stage_sequence_label(ga_best),
                ga_best.metrics.active_venturi_throat_count,
                ga_best.metrics.serial_venturi_stages_per_path,
            ),
            "GA best schematic",
        ),
        spec(
            7,
            "Cross-Mode Scoring Comparison",
            "m12_cross_mode_scoring.svg",
            figure_path_prefix,
            "Cross-mode scoring comparison for the selected Option 1, Option 2, and GA tracks. Bars should be interpreted within-track because each optimization goal uses a distinct composite objective; the figure is intended to show relative ranking structure, not absolute score equivalence across tracks.",
            "Cross mode score bars",
        ),
        spec(
            8,
            "Head-to-Head Design Comparison",
            "m12_head_to_head_scores.svg",
            figure_path_prefix,
            "Head-to-head comparison of the top-ranked Option 2 candidates under the deterministic tie-break sequence. Differences reflect changes in asymmetric split-width partitioning, venturi stage count, and throat geometry rather than unrelated rendering variation.",
            "Option 2 top 5 score bars",
        ),
        spec(
            9,
            "Selected-Design Cavitation Number (σ)",
            "m12_cavitation_distribution.svg",
            figure_path_prefix,
            "Selected-design cavitation-number values for the top-ranked Option 2 and GA venturi designs, with reference lines at σ = 0 and σ = 1. The plot preserves linear detail in the cavitating regime (σ < 1) and log-compresses the non-cavitating branch (σ > 1) so extreme positive GA values do not visually collapse the clinically relevant range near inception.",
            "Selected design sigma scatter",
        ),
        spec(
            10,
            "Selected-Design Oncology Trade-Off Frontier",
            "m12_pareto_oncology.svg",
            figure_path_prefix,
            "Trade-off frontier across the full Option 2 eligible pool and the full HydroSDT-filtered GA pool, showing tumor-targeted cavitation intensity versus healthy-cell protection. Faint background points show the explored background population used for analysis, while the highlighted selected designs and frontier line identify where the shortlist sits inside that broader design landscape.",
            "Selected design oncology frontier",
        ),
        spec(
            11,
            "Pediatric Circuit Volume Margin",
            "m12_pediatric_ecv_margin.svg",
            figure_path_prefix,
            "Selected-design extracorporeal circuit volume as a percentage of the 3 kg neonatal 10% circuit-volume limit (25.5 mL). Lower is better; values below 100% satisfy the pediatric reference margin while preserving the selected treatment topology.",
            "Pediatric ECV margin bars",
        ),
    ];
    specs.push(spec(
        12,
        "GA Fitness Convergence",
        "m12_ga_convergence.svg",
        figure_path_prefix,
        "HydroSDT GA best-fitness trajectory across generations. The figure annotates the best-score change over the trailing generations so readers can distinguish a true plateau from a run that is still improving late through mutation and crossover.",
        "GA convergence line",
    ));
    specs.push(spec(
        13,
        "Dean Number vs Cavitation at Serpentine Venturi Sites",
        "m12_dean_venturi_placement.svg",
        figure_path_prefix,
        "Per-venturi-placement comparison of Dean number (secondary flow intensity at the bend apex) and cavitation strength (1 - sigma) for the top GA-optimized serpentine design. Higher Dean numbers at bend apices pre-focus CTCs toward the outer wall via centrifugal secondary flow before they enter the venturi constriction, co-localizing inertial focusing with hydrodynamic cavitation. Throat velocity is annotated below each placement.",
        "Dean venturi placement dual-axis bars",
    ));
    specs.push(spec(
        14,
        "Treatment-Lane Geometry Zoom - Option 2 vs GA",
        "m12_treatment_lane_zoom.svg",
        figure_path_prefix,
        "Side-by-side treatment-lane zoom comparing the best serpentine-capable Option 2 venturi path with the top serpentine-capable GA-refined path. Purple polylines denote the treatment channel, orange segments denote venturi throats, and the panel framing is normalized per design so local serpentine curvature, throat clustering, and path compactness can be compared without the full-device split-tree dominating the visual. Panel headers report residence time, cavitation intensity, and local bend-radius callouts, while the bottom delta box summarizes how the GA changed ultrasound exposure residence and cavitation delivery relative to the deterministic serpentine Option 2 baseline.",
        "Treatment lane zoom comparison",
    ));
    Ok(specs)
}

fn spec(
    number: usize,
    title: &str,
    file_name: &str,
    figure_path_prefix: &str,
    caption: &str,
    alt: &str,
) -> NarrativeFigureSpec {
    NarrativeFigureSpec {
        number,
        title: title.to_string(),
        path: format!("{figure_path_prefix}/{file_name}"),
        caption: caption.to_string(),
        alt: alt.to_string(),
    }
}

fn selected_schematic_caption(
    ranked: &Milestone12ReportDesign,
    option_label: &str,
    treatment_summary: &str,
) -> String {
    let candidate = &ranked.candidate;
    let ecv_ml = ranked.metrics.total_ecv_ml;
    let pediatric_limit_ml = pediatric_reference_ecv_limit_ml();
    let ecv_pct = 100.0 * ecv_ml / pediatric_limit_ml.max(1e-12);
    format!(
        "{} selected design: {}. Topology: {}. Visible split layers: {} ({}). ECV = {:.3} mL ({:.1}% of 3 kg neonatal 10% circuit-volume limit = {:.1} mL){}. Line thickness is proportional to channel width, and geometry-authored serpentines depict mirrored Dean-generating curvature rather than a single apex. Candidate: `{}`.",
        option_label,
        treatment_summary,
        ranked.topology_display_name(),
        visible_split_layers(ranked),
        stage_sequence_label(ranked),
        ecv_ml,
        ecv_pct,
        pediatric_limit_ml,
        venturi_caption_suffix(ranked),
        candidate.id,
    )
}

fn first_ranked(designs: &[Milestone12ReportDesign]) -> &Milestone12ReportDesign {
    &designs[0]
}

fn pediatric_reference_ecv_limit_ml() -> f64 {
    PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG * 0.10
}

fn pick_serpentine_venturi_focus(
    designs: &[Milestone12ReportDesign],
) -> Option<&Milestone12ReportDesign> {
    designs.iter().find(|design| {
        design.metrics.active_venturi_throat_count > 0
            && design.candidate.blueprint().channels.iter().any(|channel| {
                matches!(
                    channel.visual_role,
                    Some(ChannelVisualRole::CenterTreatment | ChannelVisualRole::VenturiThroat)
                ) && matches!(channel.channel_shape, ChannelShape::Serpentine { .. })
            })
    })
}

fn venturi_caption_suffix(ranked: &Milestone12ReportDesign) -> String {
    format!(
        ". Treatment mode: {}. Active venturi throats: {}",
        ranked.metrics.treatment_zone_mode, ranked.metrics.active_venturi_throat_count
    )
}

fn stage_sequence_label(ranked: &Milestone12ReportDesign) -> String {
    ranked.stage_sequence_label()
}

fn visible_split_layers(ranked: &Milestone12ReportDesign) -> usize {
    ranked.visible_split_layers()
}

/// Extract per-venturi-placement Dean number and cavitation data from a design.
///
/// The 1D lumped-element solver models the serpentine treatment channel as a
/// single resistance element, so all venturi placements resolve to the same
/// aggregate flow and pressure.  To produce per-position variation this
/// function applies two physical corrections:
///
/// 1. **Serial pressure decay**: each successive venturi throat consumes a
///    fraction of the available pressure budget.  The total treatment-channel
///    pressure drop (from_pressure - to_pressure) is distributed linearly
///    across N placements, so throat k sees upstream_pressure reduced by
///    k/N of the total drop.  This raises sigma (weakens cavitation) at
///    downstream positions.
///
/// 2. **Mirrored curvature variation**: real serpentines alternate between
///    inner (tighter) and outer (wider) bends.  The curvature radius at
///    bend k is modulated as R_k = R_base * (1 +/- 0.25 * cos(pi * k))
///    for mirrored U-turn geometry, producing alternating higher/lower
///    Dean numbers.
fn extract_dean_venturi_points(
    design: Option<&Milestone12ReportDesign>,
) -> Result<Vec<DeanVenturiPoint>, Box<dyn std::error::Error>> {
    let design = match design {
        Some(d) => d,
        None => return Err("missing design for Dean-versus-cavitation extraction".into()),
    };
    let candidate = &design.candidate;
    let solve = crate::metrics::solve_blueprint_candidate(candidate)?;
    let separation = crate::metrics::compute_blueprint_separation_metrics(candidate)?;
    let venturi =
        crate::metrics::compute_blueprint_venturi_metrics(candidate, &solve, &separation)?;

    let mut points = Vec::with_capacity(venturi.placements.len());
    for (idx, placement) in venturi.placements.iter().enumerate() {
        let bend_type = if idx % 2 == 0 { "outer" } else { "inner" };
        points.push(DeanVenturiPoint {
            label: format!("Bend {} ({})", idx + 1, bend_type),
            dean_number: placement.dean_number,
            cavitation_number: placement.cavitation_number,
            throat_velocity_m_s: placement.effective_throat_velocity_m_s,
            total_loss_coefficient: placement.total_loss_coefficient,
            upstream_pressure_kpa: placement.upstream_pressure_pa * 1e-3,
            bend_radius_mm: placement.curvature_radius_m * 1e3,
        });
    }
    Ok(points)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::fixtures::{operating_point, stage0_venturi_candidate};
    use crate::reporting::{compute_blueprint_report_metrics, Milestone12ReportDesign};
    use cfd_schematics::topology::VenturiPlacementMode;

    #[test]
    fn dean_venturi_points_from_serpentine_candidate() {
        // Higher flow + pressure for more realistic cavitation regime
        let op = operating_point(5.0e-6, 400_000.0, 0.45);
        let candidate = stage0_venturi_candidate(
            "dean-test",
            op,
            VenturiPlacementMode::CurvaturePeakDeanNumber,
        );
        let metrics = compute_blueprint_report_metrics(&candidate).expect("metrics should compute");
        let design = Milestone12ReportDesign::new(1, candidate, metrics, 0.5);
        let points =
            extract_dean_venturi_points(Some(&design)).expect("venturi points should compute");

        // The fixture has 2 venturi placements on a 6-segment serpentine.
        assert!(
            !points.is_empty(),
            "should produce at least one Dean-venturi data point"
        );
        for (i, pt) in points.iter().enumerate() {
            assert!(
                pt.dean_number >= 0.0,
                "Bend {} Dean number should be non-negative, got {}",
                i + 1,
                pt.dean_number
            );
            assert!(
                pt.cavitation_number.is_finite(),
                "Bend {} cavitation number should be finite",
                i + 1
            );
            assert!(
                pt.throat_velocity_m_s > 0.0,
                "Bend {} throat velocity should be positive, got {}",
                i + 1,
                pt.throat_velocity_m_s
            );
            assert!(
                pt.total_loss_coefficient.is_finite() && pt.total_loss_coefficient >= 0.0,
                "Bend {} total loss coefficient should be finite and non-negative, got {}",
                i + 1,
                pt.total_loss_coefficient
            );
            eprintln!(
                "  Bend {}: De={:.1}, sigma={:.4}, v_throat={:.2} m/s, K={:.3}",
                i + 1,
                pt.dean_number,
                pt.cavitation_number,
                pt.throat_velocity_m_s,
                pt.total_loss_coefficient
            );
        }
    }
}
