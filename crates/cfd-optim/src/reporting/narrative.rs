//! Milestone 12 narrative report generation from template + canonical results.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::analysis::RobustnessReport;
use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::contract_trace::load_m12_contract_text;
use crate::reporting::figures::{generate_m12_report_figures, FigureGenerationInput};
use crate::reporting::narrative_addenda::{
    build_appendix_a_supplemental, build_conclusions, build_figure_sections,
    build_figure_toc_rows, build_references_block, build_storage_artifact_index,
    build_storage_policy_section,
};
use crate::reporting::narrative_sections::{
    build_cavitation_formulas_intro, build_cri_expansion_sensitivity, build_limits_of_usage,
    build_limits_of_usage_intro, build_option2_top5_table, build_results_intro,
    build_robustness_section, build_selected_table, build_selected_table_intro,
    build_strict_core_intro, build_strict_core_table, build_top5_intro, build_tri_cell_intro,
    build_tri_cell_table, build_validation_section,
};
use crate::reporting::template::render_template_strict;
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

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
    pub option1_ranked: &'a [Milestone12ReportDesign],
    pub option2_ranked: &'a [Milestone12ReportDesign],
    pub ga_top: &'a [Milestone12ReportDesign],
    pub option2_pool_all: &'a [Milestone12ReportDesign],
    pub ga_pool_all: &'a [Milestone12ReportDesign],
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

    let option1 = input.option1_ranked.first();
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
            ga_top: input.ga_top,
            option2_pool_all: input.option2_pool_all,
            ga_pool_all: input.ga_pool_all,
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
        build_storage_artifact_index(),
    );
    values.insert(
        "CONCLUSIONS".to_string(),
        build_conclusions(input.total_candidates, option1, option2, ga_best),
    );
    values.insert("REFERENCES_BLOCK".to_string(), build_references_block());

    values.insert(
        "APPENDIX_A_SUPPLEMENTAL".to_string(),
        build_appendix_a_supplemental(
            input.total_candidates,
            input.option1_pool_len,
            input.option2_pool_len,
            ga_best,
            input.option2_robustness,
            input.validation_rows,
            canonical_results_path,
        ),
    );

    let manifest_json = std::fs::read_to_string(&manifest_path)?;
    values.insert(
        "APPENDIX_B_FIGURES".to_string(),
        format!("```json\n{}\n```", manifest_json),
    );

    let narrative = render_template_strict(&template, &values)?;

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
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
) {
    let abstract_validation_sentence = if input.fast_mode {
        "Robustness screening was not computed in this fast-mode run, and multi-fidelity pressure-drop confirmation is generated by the companion `milestone12_validation` example rather than the main report run."
    } else {
        "Selected designs were screened for robustness under +/-10%/+/-20% operating-parameter perturbations. Multi-fidelity pressure-drop confirmation is generated separately by the companion `milestone12_validation` example using 1D Hagen-Poiseuille, 2D FVM, and 3D FEM solvers with non-Newtonian Carreau-Yasuda blood rheology."
    };
    let methods_pipeline_step7 = if input.fast_mode {
        "7. **Deferred robustness and companion validation flow** — When `M12_FAST=1` is enabled, robustness screening is intentionally skipped. Multi-fidelity 2D/3D validation is produced separately by the `milestone12_validation` example and may be appended after the main report run."
    } else {
        "7. **Companion multi-fidelity validation** — Selected designs are first ranked and screened here; 2D FVM and 3D FEM pressure-drop confirmation is then produced by the separate `milestone12_validation` example, while this report run retains robustness-screening results."
    };
    let milestone_scope_sentence = if input.fast_mode {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This fast-mode report documents the deterministic execution of the CFDrs pipeline from candidate generation through strict gating, scoring, and ranking, while explicitly deferring robustness screening and leaving multi-fidelity validation to the companion example."
    } else {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This report documents the deterministic execution of the CFDrs pipeline from candidate generation through strict gating, scoring, ranking, and robustness screening, with multi-fidelity validation emitted separately by the companion example."
    };
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
    values.insert(
        "OPTION1_ID".to_string(),
        option1.map_or_else(
            || "NONE_ELIGIBLE_CURRENT_PHYSICS".to_string(),
            |design| design.candidate.id.clone(),
        ),
    );
    values.insert("OPTION2_ID".to_string(), option2.candidate.id.clone());
    values.insert("GA_ID".to_string(), ga_best.candidate.id.clone());
    values.insert("GA_SEED".to_string(), M12_GA_HYDRO_SEED.to_string());
    values.insert(
        "ABSTRACT_VALIDATION_SENTENCE".to_string(),
        abstract_validation_sentence.to_string(),
    );
    values.insert(
        "METHODS_PIPELINE_STEP7".to_string(),
        methods_pipeline_step7.to_string(),
    );
    values.insert(
        "MILESTONE_SCOPE_SENTENCE".to_string(),
        milestone_scope_sentence.to_string(),
    );
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
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) {
    values.insert(
        "RESULTS_INTRO".to_string(),
        build_results_intro(
            input.total_candidates,
            input.option1_pool_len,
            input.option2_pool_len,
        ),
    );
    values.insert(
        "SELECTED_TABLE_INTRO".to_string(),
        build_selected_table_intro(option1, option2),
    );
    values.insert(
        "SELECTED_TABLE".to_string(),
        build_selected_table(option1, option2),
    );
    values.insert(
        "TOP5_INTRO".to_string(),
        build_top5_intro(input.option2_ranked),
    );
    values.insert(
        "OPTION2_TOP5_TABLE".to_string(),
        build_option2_top5_table(input.option2_ranked),
    );
    values.insert(
        "TRI_CELL_INTRO".to_string(),
        build_tri_cell_intro(option2, input.ga_top.first()),
    );
    values.insert(
        "TRI_CELL_TABLE".to_string(),
        build_tri_cell_table(input.option2_ranked, input.ga_top),
    );
    values.insert("STRICT_CORE_INTRO".to_string(), build_strict_core_intro());
    values.insert(
        "STRICT_CORE_TABLE".to_string(),
        build_strict_core_table(input.option2_ranked, input.ga_top),
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
        "CAVITATION_FORMULAS_INTRO".to_string(),
        build_cavitation_formulas_intro(option2),
    );
    values.insert(
        "LIMITS_OF_USAGE_INTRO".to_string(),
        build_limits_of_usage_intro(option2),
    );
    values.insert(
        "LIMITS_OF_USAGE".to_string(),
        build_limits_of_usage(option2),
    );
    values.insert(
        "CRI_EXPANSION_SENSITIVITY".to_string(),
        build_cri_expansion_sensitivity(),
    );
}
