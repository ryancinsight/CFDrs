//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use serde::{Deserialize, Serialize};

use super::process::{write_creation_optimization_process_figure, write_placeholder};
use super::svg::{
    write_cavitation_distribution_figure, write_cross_mode_figure, write_ga_convergence_figure,
    write_head_to_head_figure, write_pareto_figure,
    write_pediatric_ecv_figure,
};
use crate::constraints::{PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG};
use crate::reporting::{Milestone12ReportDesign, ParetoPoint};

/// Figure metadata used for dynamic table-of-contents and section rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NarrativeFigureSpec {
    pub number: usize,
    pub title: String,
    pub path: String,
    pub caption: String,
    pub alt: String,
    /// True when this figure was generated as a placeholder because the source
    /// schematic or data was not available in the current run.
    #[serde(default)]
    pub is_placeholder: bool,
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
/// # Errors
/// Returns an error if any figure cannot be written.
pub fn generate_m12_report_figures(
    figures_dir: &Path,
    figure_path_prefix: &str,
    input: &FigureGenerationInput<'_>,
) -> Result<Vec<NarrativeFigureSpec>, Box<dyn std::error::Error>> {
    std::fs::create_dir_all(figures_dir)?;

    write_creation_optimization_process_figure(
        &figures_dir.join("milestone12_creation_optimization_process.svg"),
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate.svg"),
        "Treatment Zone Layout — Bifurcation Concept",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate_trifurcation.svg"),
        "Treatment Zone Layout — Trifurcation Concept",
    )?;
    let option1_is_placeholder = if input.option1_ranked.is_empty() {
        write_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Option 1 Unavailable Under Current Physics",
            "No selective acoustic design satisfied strict eligibility under the current physics regime.",
        )?;
        true
    } else {
        ensure_existing_or_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Selected Option 1 Design — Selective Acoustic",
        )?
    };
    let option2_is_placeholder = ensure_existing_or_placeholder(
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Selected Option 2 Design — Combined Selective Venturi",
    )?;
    let ga_is_placeholder = ensure_existing_or_placeholder(
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Best HydroSDT GA-Optimized Design",
    )?;

    write_cross_mode_figure(
        &figures_dir.join("m12_cross_mode_scoring.svg"),
        input.option1_ranked,
        input.option2_ranked,
        input.ga_top,
    )?;
    write_head_to_head_figure(
        &figures_dir.join("m12_head_to_head_scores.svg"),
        input.option2_ranked,
    )?;
    write_cavitation_distribution_figure(
        &figures_dir.join("m12_cavitation_distribution.svg"),
        input.option2_ranked,
        input.ga_top,
    )?;
    write_pareto_figure(
        &figures_dir.join("m12_pareto_oncology.svg"),
        input.option2_ranked,
        input.ga_top,
        input.option2_pool_all,
        input.ga_pool_all,
    )?;
    write_pediatric_ecv_figure(
        &figures_dir.join("m12_pediatric_ecv_margin.svg"),
        input.option1_ranked.first(),
        first_ranked(input.option2_ranked),
        first_ranked(input.ga_top),
    )?;
    if !input.ga_best_per_gen.is_empty() {
        write_ga_convergence_figure(
            &figures_dir.join("m12_ga_convergence.svg"),
            input.ga_best_per_gen,
            input.fast_mode,
        )?;
    }

    let option1 = input.option1_ranked.first();
    let option2 = &input.option2_ranked[0];
    let ga_best = &input.ga_top[0];

    let make_spec = |number, title: &str, file_name, prefix, caption: &str, alt, is_ph| {
        if is_ph {
            spec_placeholder(number, title, file_name, prefix, caption, alt)
        } else {
            spec(number, title, file_name, prefix, caption, alt)
        }
    };

    let mut specs = vec![
        option1.map_or_else(
            || {
                spec_placeholder(
                    4,
                    "Option 1 Unavailable Under Current Physics",
                    "selected_option1_schematic.svg",
                    figure_path_prefix,
                    "No selective acoustic design satisfied strict eligibility under the current physics regime, so Figure 4 is an explicit placeholder documenting the empty Option 1 shortlist.",
                    "Option 1 unavailable placeholder",
                )
            },
            |option1| {
                make_spec(
                    4,
                    &format!(
                        "Selected Option 1 Design — {} Selective Acoustic",
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
                    option1_is_placeholder,
                )
            },
        ),
        make_spec(
            5,
            &format!(
                "Selected Option 2 Design — {} Combined Selective Venturi",
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
            option2_is_placeholder,
        ),
        make_spec(
            6,
            &format!(
                "Best HydroSDT GA-Optimized Design — {}",
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
            ga_is_placeholder,
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
            "Selected-design cavitation-number values for the top-ranked Option 2 and GA venturi designs, with reference lines at σ = 0 and σ = 1. This view is intended to diagnose whether each selected design is below vapor pressure, inception-capable, or above the cavitation threshold rather than merely counting regime buckets.",
            "Selected design sigma scatter",
        ),
        spec(
            10,
            "Selected-Design Oncology Trade-Off Frontier",
            "m12_pareto_oncology.svg",
            figure_path_prefix,
            "Trade-off frontier across the top-ranked Option 2 and GA designs only, showing tumor-targeted cavitation intensity versus healthy-cell protection. The frontier line highlights the nondominated selected designs instead of burying them inside a dense background pool.",
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
    if !input.ga_best_per_gen.is_empty() {
        specs.push(spec(
            12,
            "GA Fitness Convergence",
            "m12_ga_convergence.svg",
            figure_path_prefix,
            "HydroSDT GA best-fitness trajectory across generations. Flat convergence indicates the seeded parametric design was already locally optimal under the allowed in-place mutations; a rising curve indicates beneficial Dean-serpentine or throat-refinement updates were discovered.",
            "GA convergence line",
        ));
    }
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
        is_placeholder: false,
    }
}

fn spec_placeholder(
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
        caption: format!("{caption} *(Placeholder — source data not available in this run)*"),
        alt: alt.to_string(),
        is_placeholder: true,
    }
}

/// Returns `true` if a placeholder was generated (file did not already exist).
fn ensure_existing_or_placeholder(
    path: &Path,
    title: &str,
) -> Result<bool, Box<dyn std::error::Error>> {
    if path.exists() {
        Ok(false)
    } else {
        write_placeholder(
            path,
            title,
            "Source schematic not present in current run; placeholder generated.",
        )?;
        Ok(true)
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
        "{} selected design — {}. Topology: {}. Visible split layers: {} ({}). ECV = {:.3} mL ({:.1}% of 3 kg neonatal 10% circuit-volume limit = {:.1} mL){}. Line thickness is proportional to channel width, and geometry-authored serpentines depict mirrored Dean-generating curvature rather than a single apex. Candidate: `{}`.",
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
