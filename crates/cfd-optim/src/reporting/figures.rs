//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::reporting::figures_svg::{
    write_cavitation_distribution_figure, write_cross_mode_figure, write_ga_convergence_figure,
    write_head_to_head_figure, write_multifidelity_figure, write_pareto_figure, write_placeholder,
};
use crate::reporting::ValidationRow;
use crate::RankedDesign;

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
    pub option1_ranked: &'a [RankedDesign],
    pub option2_ranked: &'a [RankedDesign],
    pub rbc_ranked: &'a [RankedDesign],
    pub ga_top: &'a [RankedDesign],
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

    ensure_existing_or_placeholder(
        &figures_dir.join("milestone12_creation_optimization_process.svg"),
        "Milestone 12 Design Creation & Optimization Process",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate.svg"),
        "Treatment Zone Layout — Bifurcation Concept",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate_trifurcation.svg"),
        "Treatment Zone Layout — Trifurcation Concept",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("selected_ga_schematic.svg"),
        "Selected Option 1 Design — Ultrasound-Only",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("selected_cifx_combined_schematic.svg"),
        "Selected Option 2 Design — Unified Venturi SDT",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Best HydroSDT GA-Optimized Design",
    )?;

    write_cross_mode_figure(
        &figures_dir.join("m12_cross_mode_scoring.svg"),
        input.option1_ranked,
        input.option2_ranked,
        input.rbc_ranked,
        input.ga_top,
    )?;
    write_head_to_head_figure(
        &figures_dir.join("m12_head_to_head_scores.svg"),
        input.option2_ranked,
    )?;
    write_cavitation_distribution_figure(
        &figures_dir.join("m12_cavitation_distribution.svg"),
        input.option2_ranked,
        input.rbc_ranked,
        input.ga_top,
    )?;
    write_pareto_figure(
        &figures_dir.join("m12_pareto_oncology.svg"),
        input.option2_ranked,
        input.rbc_ranked,
        input.ga_top,
    )?;
    write_multifidelity_figure(
        &figures_dir.join("m12_multifidelity_dp.svg"),
        input.validation_rows,
        input.fast_mode,
    )?;
    write_ga_convergence_figure(
        &figures_dir.join("m12_ga_convergence.svg"),
        input.ga_best_per_gen,
        input.fast_mode,
    )?;

    Ok(vec![
        spec(
            4,
            "Selected Option 1 Design — Ultrasound-Only",
            "selected_ga_schematic.svg",
            "Selected Option 1 Design — Ultrasound-only branch network.",
            "Option 1 schematic",
        ),
        spec(
            5,
            "Selected Option 2 Design — Unified Venturi SDT",
            "selected_cifx_combined_schematic.svg",
            "Selected Option 2 Design — unified venturi SDT.",
            "Option 2 schematic",
        ),
        spec(
            6,
            "Best HydroSDT GA-Optimized Design",
            "top_hydrosdt_schematic.svg",
            "Best HydroSDT GA-optimized design.",
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
            "Multi-Fidelity Pressure-Drop Comparison",
            "m12_multifidelity_dp.svg",
            "Multi-fidelity pressure-drop comparison (1D/2D/3D).",
            "Multi fidelity pressure drop bars",
        ),
        spec(
            12,
            "GA Fitness Convergence",
            "m12_ga_convergence.svg",
            "HydroSDT GA convergence over generations.",
            "GA convergence line",
        ),
    ])
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
