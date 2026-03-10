//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::constraints::{PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG};
use crate::reporting::figures_svg::{
    write_cavitation_distribution_figure, write_creation_optimization_process_figure,
    write_cross_mode_figure, write_ga_convergence_figure, write_head_to_head_figure,
    write_multifidelity_figure, write_pareto_figure, write_pediatric_ecv_figure, write_placeholder,
};
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

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
    pub option2_pool_all: &'a [Milestone12ReportDesign],
    pub ga_pool_all: &'a [Milestone12ReportDesign],
    pub validation_rows: &'a [ValidationRow],
    pub ga_best_per_gen: &'a [f64],
    pub fast_mode: bool,
}

/// Generate report figure assets and return ordered figure specs.
///
/// # Errors
/// Returns an error if any figure cannot be written.
pub fn generate_m12_report_figures(
    figures_dir: &Path,
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
    if input.option1_ranked.is_empty() {
        write_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Option 1 Unavailable Under Current Physics",
            "No selective acoustic design satisfied strict eligibility under the current physics regime.",
        )?;
    } else {
        ensure_existing_or_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Selected Option 1 Design — Selective Acoustic",
        )?;
    }
    ensure_existing_or_placeholder(
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Selected Option 2 Design — Combined Selective Venturi",
    )?;
    ensure_existing_or_placeholder(
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
    if !input.validation_rows.is_empty() {
        write_multifidelity_figure(
            &figures_dir.join("m12_multifidelity_dp.svg"),
            input.validation_rows,
            input.fast_mode,
        )?;
    }
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

    let mut specs = vec![
        option1.map_or_else(
            || {
                spec(
                    4,
                    "Option 1 Unavailable Under Current Physics",
                    "selected_option1_schematic.svg",
                    "No selective acoustic design satisfied strict eligibility under the current physics regime, so Figure 4 is an explicit placeholder documenting the empty Option 1 shortlist.",
                    "Option 1 unavailable placeholder",
                )
            },
            |option1| {
                spec(
                    4,
                    &format!(
                        "Selected Option 1 Design — {} Selective Acoustic",
                        stage_sequence_label(option1)
                    ),
                    "selected_option1_schematic.svg",
                    &selected_schematic_caption(
                        option1,
                        "Option 1",
                        "selective acoustic center-treatment network",
                    ),
                    "Option 1 schematic",
                )
            },
        ),
        spec(
            5,
            &format!(
                "Selected Option 2 Design — {} Combined Selective Venturi",
                stage_sequence_label(option2)
            ),
            "selected_option2_combined_schematic.svg",
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
                "Best HydroSDT GA-Optimized Design — {}",
                stage_sequence_label(ga_best)
            ),
            "top_hydrosdt_schematic.svg",
            &format!(
                "Best HydroSDT GA-optimized design. Topology: {}. Visible split layers: {}. Active venturi throats: {}.",
                ga_best.topology_display_name(),
                visible_split_layers(ga_best),
                ga_best.metrics.active_venturi_throat_count,
            ),
            "GA best schematic",
        ),
        spec(
            7,
            "Cross-Mode Scoring Comparison",
            "m12_cross_mode_scoring.svg",
            "Cross-mode scoring comparison for selected tracks.",
            "Cross mode score bars",
        ),
        spec(
            8,
            "Head-to-Head Design Comparison",
            "m12_head_to_head_scores.svg",
            "Head-to-head score comparison for Option 2 top-ranked designs.",
            "Option 2 top 5 score bars",
        ),
        spec(
            9,
            "Cavitation Number (σ) Distribution",
            "m12_cavitation_distribution.svg",
            "Cavitation number category distribution across selected venturi designs.",
            "Cavitation sigma category bars",
        ),
        spec(
            10,
            "Pareto Front — Oncology Objectives",
            "m12_pareto_oncology.svg",
            "Pareto view of cancer-targeted cavitation versus therapeutic window score.",
            "Pareto oncology scatter",
        ),
        spec(
            11,
            "Pediatric Circuit Volume Margin",
            "m12_pediatric_ecv_margin.svg",
            "Selected-design ECV as a percentage of the 3 kg neonatal 10% circuit-volume limit (25.5 mL). Lower is better; values below 100% satisfy the pediatric reference margin.",
            "Pediatric ECV margin bars",
        ),
    ];
    if !input.validation_rows.is_empty() {
        specs.push(spec(
            12,
            "Multi-Fidelity Pressure-Drop Comparison",
            "m12_multifidelity_dp.svg",
            "Multi-fidelity pressure-drop comparison (1D/2D/3D).",
            "Multi fidelity pressure drop bars",
        ));
    }
    if !input.ga_best_per_gen.is_empty() {
        // Renumber dynamically to follow the pediatric ECV figure.
        let ga_fig_num = if input.validation_rows.is_empty() {
            12
        } else {
            13
        };
        specs.push(spec(
            ga_fig_num,
            "GA Fitness Convergence",
            "m12_ga_convergence.svg",
            "HydroSDT GA convergence over generations.",
            "GA convergence line",
        ));
    }
    Ok(specs)
}

fn spec(
    number: usize,
    title: &str,
    file_name: &str,
    caption: &str,
    alt: &str,
) -> NarrativeFigureSpec {
    NarrativeFigureSpec {
        number,
        title: title.to_string(),
        path: format!("../report/figures/{file_name}"),
        caption: caption.to_string(),
        alt: alt.to_string(),
    }
}

fn ensure_existing_or_placeholder(
    path: &Path,
    title: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if path.exists() {
        Ok(())
    } else {
        write_placeholder(
            path,
            title,
            "Source schematic not present in current run; placeholder generated.",
        )
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
        "{} selected design — {}. Topology: {}. Visible split layers: {} ({}). ECV = {:.3} mL ({:.1}% of 3 kg neonatal 10% circuit-volume limit = {:.1} mL){}. Candidate: `{}`.",
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
