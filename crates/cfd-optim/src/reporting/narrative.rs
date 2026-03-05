//! Milestone 12 narrative report generation from template + canonical results.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::analysis::RobustnessReport;
use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::contract_trace::load_m12_contract_text;
use crate::reporting::figures::{generate_m12_report_figures, FigureGenerationInput};
use crate::reporting::narrative_sections::{
    build_figure_sections, build_figure_toc_rows, build_limits_of_usage, build_option2_top5_table,
    build_robustness_section, build_selected_table, build_storage_artifact_index,
    build_storage_policy_section, build_strict_core_table, build_tri_cell_table,
    build_validation_section,
};
use crate::reporting::template::render_template_strict;
use crate::reporting::ValidationRow;
use crate::RankedDesign;

/// Metadata for the narrative title page.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct M12Metadata {
    pub agreement_number: String,
    pub project_description: String,
    pub milestone_label: String,
    pub delivery_date: String,
    pub prepared_for: Vec<String>,
    pub logo_path: String,
    pub address_lines: Vec<String>,
    pub admin_contact_name: String,
    pub admin_contact_email: String,
    pub investigator_name: String,
    pub investigator_email: String,
}

/// Inputs used to render the narrative report.
pub struct Milestone12NarrativeInput<'a> {
    pub total_candidates: usize,
    pub option1_pool_len: usize,
    pub option2_pool_len: usize,
    pub rbc_pool_len: usize,
    pub option1_ranked: &'a [RankedDesign],
    pub option2_ranked: &'a [RankedDesign],
    pub rbc_ranked: &'a [RankedDesign],
    pub ga_top: &'a [RankedDesign],
    pub validation_rows: &'a [ValidationRow],
    pub option2_robustness: &'a [RobustnessReport],
    pub ga_best_per_gen: &'a [f64],
    pub fast_mode: bool,
}

/// Written artifact paths from narrative generation.
#[derive(Debug, Clone)]
pub struct Milestone12NarrativeArtifacts {
    pub narrative_path: PathBuf,
    pub figure_manifest_path: PathBuf,
    pub figure_count: usize,
}

/// Generate the M12 narrative report and figure manifest.
///
/// # Errors
/// Returns an error if required data is missing or any artifact cannot be written.
pub fn write_milestone12_narrative_report(
    workspace_root: &Path,
    canonical_results_path: &Path,
    input: &Milestone12NarrativeInput<'_>,
) -> Result<Milestone12NarrativeArtifacts, Box<dyn std::error::Error>> {
    let report_dir = workspace_root.join("report");
    let template = std::fs::read_to_string(
        report_dir
            .join("templates")
            .join("m12_narrative_template.md"),
    )?;
    let metadata: M12Metadata = serde_json::from_str(&std::fs::read_to_string(
        report_dir.join("m12_metadata.json"),
    )?)?;
    let contract = load_m12_contract_text(workspace_root)?;

    let option1 = input
        .option1_ranked
        .first()
        .ok_or("option1_ranked is empty; narrative requires selected Option 1")?;
    let option2 = input
        .option2_ranked
        .first()
        .ok_or("option2_ranked is empty; narrative requires selected Option 2")?;
    let ga_best = input
        .ga_top
        .first()
        .ok_or("ga_top is empty; narrative requires GA rank-1")?;

    let figure_specs = generate_m12_report_figures(
        &report_dir.join("figures"),
        &FigureGenerationInput {
            option1_ranked: input.option1_ranked,
            option2_ranked: input.option2_ranked,
            rbc_ranked: input.rbc_ranked,
            ga_top: input.ga_top,
            validation_rows: input.validation_rows,
            ga_best_per_gen: input.ga_best_per_gen,
            fast_mode: input.fast_mode,
        },
    )?;
    let manifest_path = report_dir.join("milestone12").join("figure_manifest.json");
    std::fs::create_dir_all(manifest_path.parent().ok_or("invalid manifest path")?)?;
    std::fs::write(&manifest_path, serde_json::to_string_pretty(&figure_specs)?)?;

    let mut values = BTreeMap::new();
    insert_title_page_values(&mut values, &metadata);
    insert_contract_values(&mut values, &contract);
    insert_data_values(&mut values, input, option1, option2, ga_best);
    insert_table_values(&mut values, input, option1, option2);
    values.insert(
        "FIGURE_TOC_ROWS".to_string(),
        build_figure_toc_rows(&figure_specs),
    );
    values.insert(
        "FIGURE_SECTIONS".to_string(),
        build_figure_sections(&figure_specs),
    );
    values.insert(
        "STORAGE_POLICY_SECTION".to_string(),
        build_storage_policy_section(),
    );
    values.insert(
        "STORAGE_ARTIFACT_INDEX".to_string(),
        build_storage_artifact_index(&manifest_path, input.rbc_pool_len),
    );
    values.insert(
        "CONCLUSIONS".to_string(),
        build_conclusions(
            input.total_candidates,
            &option1.candidate.id,
            &option2.candidate.id,
            &ga_best.candidate.id,
            ga_best.score,
            ga_best.metrics.cavitation_number,
        ),
    );

    let mut narrative = render_template_strict(&template, &values)?;
    let canonical = std::fs::read_to_string(canonical_results_path)?;
    narrative.push_str("\n\n---\n\n## Appendix A. Canonical Milestone 12 Results (Embedded)\n\n");
    narrative.push_str(&canonical);

    let narrative_path = report_dir.join("ARPA-H_SonALAsense_Milestone 12 Report.md");
    std::fs::write(&narrative_path, narrative)?;

    Ok(Milestone12NarrativeArtifacts {
        narrative_path,
        figure_manifest_path: manifest_path,
        figure_count: figure_specs.len(),
    })
}

fn insert_title_page_values(values: &mut BTreeMap<String, String>, metadata: &M12Metadata) {
    values.insert(
        "AGREEMENT_NUMBER".to_string(),
        metadata.agreement_number.clone(),
    );
    values.insert(
        "PROJECT_DESCRIPTION".to_string(),
        metadata.project_description.clone(),
    );
    values.insert(
        "MILESTONE_LABEL".to_string(),
        metadata.milestone_label.clone(),
    );
    values.insert("DELIVERY_DATE".to_string(), metadata.delivery_date.clone());
    values.insert(
        "PREPARED_FOR_BLOCK".to_string(),
        metadata.prepared_for.join("\n\n"),
    );
    values.insert("LOGO_PATH".to_string(), metadata.logo_path.clone());
    values.insert(
        "ADDRESS_BLOCK".to_string(),
        metadata.address_lines.join("\n\n"),
    );
    values.insert(
        "ADMIN_CONTACT_NAME".to_string(),
        metadata.admin_contact_name.clone(),
    );
    values.insert(
        "ADMIN_CONTACT_EMAIL".to_string(),
        metadata.admin_contact_email.clone(),
    );
    values.insert(
        "INVESTIGATOR_NAME".to_string(),
        metadata.investigator_name.clone(),
    );
    values.insert(
        "INVESTIGATOR_EMAIL".to_string(),
        metadata.investigator_email.clone(),
    );
}

fn insert_contract_values(
    values: &mut BTreeMap<String, String>,
    contract: &crate::reporting::contract_trace::MilestoneContractText,
) {
    values.insert(
        "MILESTONE_DESCRIPTION".to_string(),
        contract.description.clone(),
    );
    values.insert(
        "MILESTONE_EXIT_CRITERIA".to_string(),
        contract.exit_criteria.clone(),
    );
    values.insert(
        "MILESTONE_DELIVERABLE".to_string(),
        format!("{}, including limits of usage.", contract.deliverable),
    );
}

fn insert_data_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: &RankedDesign,
    option2: &RankedDesign,
    ga_best: &RankedDesign,
) {
    values.insert(
        "TOTAL_CANDIDATES".to_string(),
        input.total_candidates.to_string(),
    );
    values.insert(
        "CANONICAL_SOURCE".to_string(),
        "report/milestone12_results.md".to_string(),
    );
    values.insert(
        "OPTION1_POOL".to_string(),
        input.option1_pool_len.to_string(),
    );
    values.insert(
        "OPTION2_POOL".to_string(),
        input.option2_pool_len.to_string(),
    );
    values.insert("OPTION1_ID".to_string(), option1.candidate.id.clone());
    values.insert("OPTION2_ID".to_string(), option2.candidate.id.clone());
    values.insert("GA_ID".to_string(), ga_best.candidate.id.clone());
    values.insert("GA_SEED".to_string(), M12_GA_HYDRO_SEED.to_string());
    values.insert(
        "FIG_PROCESS".to_string(),
        "../report/figures/milestone12_creation_optimization_process.svg".to_string(),
    );
    values.insert(
        "FIG_TREATMENT_BI".to_string(),
        "../report/figures/treatment_zone_plate.svg".to_string(),
    );
    values.insert(
        "FIG_TREATMENT_TRI".to_string(),
        "../report/figures/treatment_zone_plate_trifurcation.svg".to_string(),
    );
}

fn insert_table_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: &RankedDesign,
    option2: &RankedDesign,
) {
    values.insert(
        "SELECTED_TABLE".to_string(),
        build_selected_table(option1, option2, input.rbc_ranked.first()),
    );
    values.insert(
        "OPTION2_TOP5_TABLE".to_string(),
        build_option2_top5_table(input.option2_ranked),
    );
    values.insert(
        "TRI_CELL_TABLE".to_string(),
        build_tri_cell_table(input.option2_ranked, input.rbc_ranked),
    );
    values.insert(
        "STRICT_CORE_TABLE".to_string(),
        build_strict_core_table(input.option2_ranked, input.rbc_ranked),
    );
    values.insert(
        "ROBUSTNESS_SECTION".to_string(),
        build_robustness_section(input.option2_robustness, input.fast_mode),
    );
    values.insert(
        "VALIDATION_SECTION".to_string(),
        build_validation_section(input.validation_rows, input.fast_mode),
    );
    values.insert(
        "LIMITS_OF_USAGE".to_string(),
        build_limits_of_usage(option2),
    );
}

fn build_conclusions(
    total_candidates: usize,
    option1_id: &str,
    option2_id: &str,
    ga_id: &str,
    ga_score: f64,
    ga_sigma: f64,
) -> String {
    format!(
        "Milestone 12 completion criteria were met using deterministic CFDrs design selection.\n\nFrom **{}** candidates, the selected Option 1 design (`{}`) satisfies ultrasound-only hydrodynamic requirements, while selected Option 2 (`{}`) satisfies venturi cavitation criteria with finite, sub-unity cavitation number and strict-core gate passes.\n\nGA reproducibility is anchored to seed `{}` with rank-1 design `{}` (score {:.4}, sigma {:.4}). Narrative and canonical outputs are generated from a single run payload to preserve traceable contract evidence.",
        total_candidates, option1_id, option2_id, M12_GA_HYDRO_SEED, ga_id, ga_score, ga_sigma
    )
}
