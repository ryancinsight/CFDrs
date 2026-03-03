//! Milestone 12 two-concept selection pipeline.
//!
//! Generates report-ready ranked outputs for:
//! 1. Ultrasound branch-network concept (no venturi throats).
//! 2. Venturi CIF concept (hydrodynamic cavitation with RBC protection).
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example m12_two_concept_selection --no-default-features
//! ```

use cfd_optim::{
    save_schematic_svg, save_top5_json, OptimMode, RankedDesign, SdtOptimizer, SdtWeights,
};
use serde::Serialize;
use std::fmt::Write as _;
use std::path::Path;

#[derive(Debug, Serialize)]
struct ConceptTopline<'a> {
    concept: &'a str,
    mode: &'a str,
    candidate_id: &'a str,
    score: f64,
    topology: String,
    cavitation_number: f64,
    rbc_venturi_exposure_fraction: f64,
    cancer_center_fraction: f64,
    therapeutic_window_score: f64,
    hemolysis_index_per_pass: f64,
    wall_shear_p95_pa: f64,
    mean_residence_time_s: f64,
    total_ecv_ml: f64,
    therapy_channel_fraction: f64,
}

fn require_at_least(
    ranked: Vec<RankedDesign>,
    n: usize,
    label: &str,
) -> Result<Vec<RankedDesign>, Box<dyn std::error::Error>> {
    if ranked.len() < n {
        return Err(format!(
            "{label}: required at least {n} candidates, found {}",
            ranked.len()
        )
        .into());
    }
    Ok(ranked.into_iter().take(n).collect())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let out_dir = Path::new("report/milestone12");
    std::fs::create_dir_all(out_dir)?;

    let weights = SdtWeights::default();

    // Concept 1: ultrasound-only branch-network (no venturi).
    let exposure_ranked = SdtOptimizer::new(OptimMode::UniformExposure, weights).all_ranked()?;
    let option1_ranked = require_at_least(
        exposure_ranked
            .into_iter()
            .filter(|d| !d.candidate.topology.has_venturi())
            .collect(),
        5,
        "option1_ultrasound_branch",
    )?;

    // Concept 2: venturi-CIF with combined SDT + leukapheresis gating.
    let combined_mode = OptimMode::CombinedSdtLeukapheresis {
        leuka_weight: 0.5,
        sdt_weight: 0.5,
        patient_weight_kg: 3.0,
    };
    let venturi_ranked = SdtOptimizer::new(combined_mode, weights).all_ranked()?;
    let option2_ranked = require_at_least(
        venturi_ranked
            .into_iter()
            .filter(|d| d.candidate.topology.has_venturi())
            .collect(),
        5,
        "option2_venturi_cif",
    )?;

    save_top5_json(
        &option1_ranked,
        &out_dir.join("two_concept_option1_ultrasound_top5.json"),
    )?;
    save_top5_json(
        &option2_ranked,
        &out_dir.join("two_concept_option2_venturi_top5.json"),
    )?;

    save_schematic_svg(
        &option1_ranked[0].candidate,
        &out_dir.join("two_concept_option1_selected.svg"),
    )?;
    save_schematic_svg(
        &option2_ranked[0].candidate,
        &out_dir.join("two_concept_option2_selected.svg"),
    )?;

    let toplines = vec![
        ConceptTopline {
            concept: "Option 1: Ultrasound branch-network (no venturi)",
            mode: "UniformExposure",
            candidate_id: &option1_ranked[0].candidate.id,
            score: option1_ranked[0].score,
            topology: option1_ranked[0].candidate.topology.short().to_string(),
            cavitation_number: option1_ranked[0].metrics.cavitation_number,
            rbc_venturi_exposure_fraction: option1_ranked[0].metrics.rbc_venturi_exposure_fraction,
            cancer_center_fraction: option1_ranked[0].metrics.cancer_center_fraction,
            therapeutic_window_score: option1_ranked[0].metrics.therapeutic_window_score,
            hemolysis_index_per_pass: option1_ranked[0].metrics.hemolysis_index_per_pass,
            wall_shear_p95_pa: option1_ranked[0].metrics.wall_shear_p95_pa,
            mean_residence_time_s: option1_ranked[0].metrics.mean_residence_time_s,
            total_ecv_ml: option1_ranked[0].metrics.total_ecv_ml,
            therapy_channel_fraction: option1_ranked[0].metrics.therapy_channel_fraction,
        },
        ConceptTopline {
            concept: "Option 2: Venturi CIF/RBC-protected SDT",
            mode: "CombinedSdtLeukapheresis",
            candidate_id: &option2_ranked[0].candidate.id,
            score: option2_ranked[0].score,
            topology: option2_ranked[0].candidate.topology.short().to_string(),
            cavitation_number: option2_ranked[0].metrics.cavitation_number,
            rbc_venturi_exposure_fraction: option2_ranked[0].metrics.rbc_venturi_exposure_fraction,
            cancer_center_fraction: option2_ranked[0].metrics.cancer_center_fraction,
            therapeutic_window_score: option2_ranked[0].metrics.therapeutic_window_score,
            hemolysis_index_per_pass: option2_ranked[0].metrics.hemolysis_index_per_pass,
            wall_shear_p95_pa: option2_ranked[0].metrics.wall_shear_p95_pa,
            mean_residence_time_s: option2_ranked[0].metrics.mean_residence_time_s,
            total_ecv_ml: option2_ranked[0].metrics.total_ecv_ml,
            therapy_channel_fraction: option2_ranked[0].metrics.therapy_channel_fraction,
        },
    ];

    let topline_json = serde_json::to_string_pretty(&toplines)?;
    std::fs::write(
        out_dir.join("two_concept_selection_summary.json"),
        topline_json,
    )?;

    let mut md = String::new();
    md.push_str("# Milestone 12 Two-Concept Selection\n\n");
    md.push_str(
        "This run enforces the two Milestone 12 concept tracks on the shared \
         96-well footprint (127.76 x 85.47 mm):\n\n",
    );
    md.push_str("- Option 1: branch-network treatment channels without venturi throats (acoustic-only treatment zone).\n");
    md.push_str("- Option 2: venturi-enabled CIF treatment channels with RBC-protective peripheralization.\n\n");
    md.push_str("Top-5 JSON and selected schematics are exported to `report/milestone12/`.\n\n");
    md.push_str(
        "| Concept | Mode | Candidate | Score | Topology | Sigma | RBC venturi exposure | Cancer center frac | Therapeutic window | HI/pass | P95 shear (Pa) | Residence (s) | ECV (mL) |\n",
    );
    md.push_str(
        "| --- | --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |\n",
    );
    for row in &toplines {
        let _ = writeln!(
            md,
            "| {} | {} | `{}` | {:.4} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {:.2e} | {:.2} | {:.4} | {:.4} |",
            row.concept,
            row.mode,
            row.candidate_id,
            row.score,
            row.topology,
            row.cavitation_number,
            row.rbc_venturi_exposure_fraction,
            row.cancer_center_fraction,
            row.therapeutic_window_score,
            row.hemolysis_index_per_pass,
            row.wall_shear_p95_pa,
            row.mean_residence_time_s,
            row.total_ecv_ml,
        );
    }

    std::fs::write(out_dir.join("two_concept_selection.md"), md)?;

    println!("Saved:");
    println!("  - two_concept_option1_ultrasound_top5.json");
    println!("  - two_concept_option2_venturi_top5.json");
    println!("  - two_concept_option1_selected.svg");
    println!("  - two_concept_option2_selected.svg");
    println!("  - two_concept_selection_summary.json");
    println!("  - two_concept_selection.md");

    Ok(())
}
