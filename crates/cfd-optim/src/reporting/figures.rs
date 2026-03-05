//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::reporting::figures_svg::{
    write_cavitation_distribution_figure, write_cross_mode_figure, write_ga_convergence_figure,
    write_head_to_head_figure, write_multifidelity_figure, write_pareto_figure, write_placeholder,
};
use crate::reporting::ValidationRow;
use crate::RankedDesign;
use crate::{DesignCandidate, DesignTopology, TreatmentZoneMode};

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
        "Selected Option 1 Design — Selective Acoustic",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("selected_cifx_combined_schematic.svg"),
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

    let option1 = &input.option1_ranked[0];
    let option2 = &input.option2_ranked[0];
    let ga_best = &input.ga_top[0];

    let mut specs = vec![
        spec(
            4,
            &format!(
                "Selected Option 1 Design — {} Selective Acoustic",
                stage_sequence_label(&option1.candidate)
            ),
            "selected_ga_schematic.svg",
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
                "Selected Option 2 Design — {} Combined Selective Venturi",
                stage_sequence_label(&option2.candidate)
            ),
            "selected_cifx_combined_schematic.svg",
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
                stage_sequence_label(&ga_best.candidate)
            ),
            "top_hydrosdt_schematic.svg",
            &format!(
                "Best HydroSDT GA-optimized design. Topology: {}. Visible split layers: {}. Active venturi throats: {}.",
                ga_best.candidate.topology.name(),
                visible_split_layers(&ga_best.candidate),
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
    ];
    if !input.validation_rows.is_empty() {
        specs.push(spec(
            11,
            "Multi-Fidelity Pressure-Drop Comparison",
            "m12_multifidelity_dp.svg",
            "Multi-fidelity pressure-drop comparison (1D/2D/3D).",
            "Multi fidelity pressure drop bars",
        ));
    }
    if !input.ga_best_per_gen.is_empty() {
        specs.push(spec(
            12,
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
    ranked: &RankedDesign,
    option_label: &str,
    treatment_summary: &str,
) -> String {
    let candidate = &ranked.candidate;
    format!(
        "{} selected design — {}. Topology: {}. Visible split layers: {} ({}){}. Candidate: `{}`.",
        option_label,
        treatment_summary,
        candidate.topology.name(),
        visible_split_layers(candidate),
        stage_sequence_label(candidate),
        venturi_caption_suffix(ranked),
        candidate.id,
    )
}

fn venturi_caption_suffix(ranked: &RankedDesign) -> String {
    let treatment_mode = ranked.candidate.treatment_zone_mode_effective();
    match treatment_mode {
        TreatmentZoneMode::UltrasoundOnly => {
            format!(
                ". Treatment mode: UltrasoundOnly. Active venturi throats: {}",
                ranked.metrics.active_venturi_throat_count
            )
        }
        TreatmentZoneMode::VenturiThroats => {
            format!(
                ". Treatment mode: VenturiThroats. Active venturi throats: {}",
                ranked.metrics.active_venturi_throat_count
            )
        }
    }
}

fn stage_sequence_label(candidate: &DesignCandidate) -> String {
    match candidate.topology {
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
            vec!["Tri"; usize::from(n_levels.min(4))].join("→")
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            let mut stages = vec!["Tri"; usize::from(n_pretri.min(2))];
            stages.push("Tri");
            stages.push("Bi");
            stages.join("→")
        }
        DesignTopology::TriBiTriSelectiveVenturi => "Tri→Bi→Tri".to_owned(),
        DesignTopology::TrifurcationBifurcationVenturi => "Tri→Bi".to_owned(),
        DesignTopology::BifurcationTrifurcationVenturi => "Bi→Tri".to_owned(),
        DesignTopology::TrifurcationBifurcationBifurcationVenturi => "Tri→Bi→Bi".to_owned(),
        DesignTopology::TripleTrifurcationVenturi => "Tri→Tri→Tri".to_owned(),
        DesignTopology::QuadTrifurcationVenturi => "Tri→Tri→Tri→Tri".to_owned(),
        DesignTopology::DoubleTrifurcationVenturi => "Tri→Tri".to_owned(),
        DesignTopology::DoubleBifurcationVenturi => "Bi→Bi".to_owned(),
        DesignTopology::TripleBifurcationVenturi => "Bi→Bi→Bi".to_owned(),
        DesignTopology::BifurcationVenturi | DesignTopology::BifurcationSerpentine => {
            "Bi".to_owned()
        }
        DesignTopology::TrifurcationVenturi | DesignTopology::TrifurcationSerpentine => {
            "Tri".to_owned()
        }
        _ => candidate.topology.short().to_owned(),
    }
}

fn visible_split_layers(candidate: &DesignCandidate) -> usize {
    match candidate.topology {
        DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => {
            usize::from(n_levels.min(4))
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => {
            usize::from(n_pretri.min(2)) + 2
        }
        DesignTopology::TriBiTriSelectiveVenturi => 3,
        DesignTopology::TrifurcationBifurcationVenturi
        | DesignTopology::BifurcationTrifurcationVenturi
        | DesignTopology::DoubleBifurcationVenturi
        | DesignTopology::DoubleTrifurcationVenturi
        | DesignTopology::DoubleBifurcationSerpentine => 2,
        DesignTopology::TripleTrifurcationVenturi
        | DesignTopology::TripleBifurcationVenturi
        | DesignTopology::TrifurcationBifurcationBifurcationVenturi => 3,
        DesignTopology::QuadTrifurcationVenturi => 4,
        DesignTopology::BifurcationVenturi
        | DesignTopology::TrifurcationVenturi
        | DesignTopology::BifurcationSerpentine
        | DesignTopology::TrifurcationSerpentine
        | DesignTopology::AsymmetricBifurcationSerpentine
        | DesignTopology::AsymmetricTrifurcationVenturi
        | DesignTopology::CellSeparationVenturi
        | DesignTopology::WbcCancerSeparationVenturi => 1,
        _ => 0,
    }
}
