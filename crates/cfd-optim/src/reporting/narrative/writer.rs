//! Milestone 12 narrative report generation from template + canonical results.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::super::asset_review::write_asset_review_manifest;
use super::super::validation_traceability::{
    build_validation_summary_table, build_validation_traceability_intro,
};
use super::addenda::{
    build_conclusions, build_figure_sections, build_figure_toc_rows, build_references_block,
    build_storage_artifact_index, build_storage_policy_section,
    build_workspace_configuration_section,
};
use super::contract::load_m12_contract_text;
use super::sections::{
    build_cavitation_formulas_intro, build_cri_expansion_sensitivity, build_limits_of_usage,
    build_limits_of_usage_intro, build_option2_top5_table, build_results_intro,
    build_robustness_section, build_selected_table, build_selected_table_intro,
    build_serpentine_venturi_discussion, build_strict_core_intro, build_strict_core_table,
    build_top5_intro, build_treatment_time_analysis, build_tri_cell_intro, build_tri_cell_table,
    physics_model_inventory,
};
use super::template::render_template_strict;
use crate::analysis::RobustnessReport;
use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::figures::{generate_m12_report_figures, FigureGenerationInput};
use crate::reporting::{Milestone12ReportDesign, ParetoPoint, ValidationRow};

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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Milestone12GaRankingAuditEntry {
    pub rank: usize,
    pub candidate_id: String,
    pub score: f64,
    pub adjusted_selection_score: f64,
    pub geometry_concentration_penalty: f64,
    pub operating_point_diversity_penalty: f64,
}

/// Inputs used to render the narrative report.
pub struct Milestone12NarrativeInput<'a> {
    pub authoritative_run: bool,
    pub canonical_source: String,
    pub total_candidates: usize,
    pub option1_evaluated_count: usize,
    pub option2_evaluated_count: usize,
    pub option1_pool_len: usize,
    pub option2_pool_len: usize,
    pub option1_sequence_summary_markdown: String,
    pub option1_ranked: &'a [Milestone12ReportDesign],
    pub option2_ranked: &'a [Milestone12ReportDesign],
    pub ga_top: &'a [Milestone12ReportDesign],
    pub option2_pool_all: &'a [ParetoPoint],
    pub ga_pool_all: &'a [ParetoPoint],
    pub validation_rows: &'a [ValidationRow],
    pub option2_robustness: &'a [RobustnessReport],
    pub ga_best_per_gen: &'a [f64],
    pub ga_ranking_audit: &'a [Milestone12GaRankingAuditEntry],
    pub topology_family_count: usize,
    pub fast_mode: bool,
}

/// Written artifact paths from narrative generation.
#[derive(Debug, Clone)]
pub struct Milestone12NarrativeArtifacts {
    pub narrative_path: PathBuf,
    pub figure_manifest_path: PathBuf,
    pub asset_review_manifest_path: PathBuf,
    pub asset_review_complete: bool,
    pub figure_count: usize,
}

/// Generate the M12 narrative report and figure manifest.
///
/// # Errors
/// Returns an error if required data is missing or any artifact cannot be written.
pub fn write_milestone12_narrative_report(
    workspace_root: &Path,
    _canonical_results_path: &Path,
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

    let artifact_root = report_dir.clone();
    let figures_dir = report_dir.join("figures");
    let figure_path_prefix = "../report/figures";

    let figure_specs = generate_m12_report_figures(
        &figures_dir,
        figure_path_prefix,
        &FigureGenerationInput {
            option1_ranked: input.option1_ranked,
            option2_ranked: input.option2_ranked,
            ga_top: input.ga_top,
            option2_pool_all: input.option2_pool_all,
            ga_pool_all: input.ga_pool_all,
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
    let asset_review_manifest_path = report_dir
        .join("milestone12")
        .join("asset_review_manifest.json");
    let validation_rows_path = report_dir
        .join("milestone12")
        .join("milestone12_validation_rows.json");
    insert_data_values(
        &mut values,
        input,
        option1,
        option2,
        ga_best,
        figure_specs.len(),
        &manifest_path,
        &asset_review_manifest_path,
        &validation_rows_path,
    );
    insert_table_values(&mut values, input, option1, option2);
    let narrative_only_specs = figure_specs
        .iter()
        .filter(|spec| spec.number >= 4)
        .cloned()
        .collect::<Vec<_>>();
    values.insert(
        "FIGURE_TOC_ROWS".to_string(),
        build_figure_toc_rows(&narrative_only_specs),
    );
    values.insert(
        "FIGURE_SECTIONS".to_string(),
        build_figure_sections(&narrative_only_specs),
    );
    let figure_count = figure_specs.len();
    drop(narrative_only_specs);
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
        build_conclusions(
            input.total_candidates,
            input.topology_family_count,
            option1,
            option2,
            input.ga_top,
            ga_best,
            input.ga_best_per_gen,
            input.ga_ranking_audit,
            if input.option2_robustness.is_empty() {
                16
            } else {
                17
            },
        ),
    );
    values.insert(
        "WORKSPACE_CONFIGURATION_SECTION".to_string(),
        build_workspace_configuration_section(),
    );
    values.insert(
        "TABLE_DERIVED_CAPTION_NUMBER".to_string(),
        if input.option2_robustness.is_empty() {
            "14".to_string()
        } else {
            "15".to_string()
        },
    );
    values.insert(
        "TABLE_ABBREVIATIONS_CAPTION_NUMBER".to_string(),
        if input.option2_robustness.is_empty() {
            "17".to_string()
        } else {
            "18".to_string()
        },
    );
    values.insert("REFERENCES_BLOCK".to_string(), build_references_block());

    let narrative = normalize_generated_markdown(&render_template_strict(&template, &values)?);
    // Template string and BTreeMap values consumed — free ~300 KB.
    drop(template);
    drop(values);

    std::fs::create_dir_all(&artifact_root)?;
    let narrative_path = report_dir.join("ARPA-H_SonALAsense_Milestone 12 Report.md");
    std::fs::write(&narrative_path, narrative)?;
    let asset_review = write_asset_review_manifest(&report_dir, &narrative_path, &figure_specs)?;
    drop(figure_specs);

    Ok(Milestone12NarrativeArtifacts {
        narrative_path,
        figure_manifest_path: manifest_path,
        asset_review_manifest_path: asset_review.manifest_path,
        asset_review_complete: asset_review.complete,
        figure_count,
    })
}

fn normalize_generated_markdown(input: &str) -> String {
    let mut normalized = String::with_capacity(input.len());
    let mut blank_run = 0usize;

    for line in input.lines() {
        if line.trim().is_empty() {
            blank_run += 1;
            if blank_run > 1 {
                continue;
            }
            normalized.push('\n');
            continue;
        }

        blank_run = 0;
        normalized.push_str(line);
        normalized.push('\n');
    }

    while normalized.ends_with("\n\n") {
        normalized.pop();
    }

    if !input.ends_with('\n') && !normalized.is_empty() {
        normalized.pop();
    } else if input.ends_with('\n') && !normalized.ends_with('\n') {
        normalized.push('\n');
    }

    normalized
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
    contract: &super::contract::MilestoneContractText,
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
        format_milestone_deliverable(&contract.deliverable),
    );
}

fn format_milestone_deliverable(deliverable: &str) -> String {
    let trimmed = deliverable.trim();
    if trimmed.to_ascii_lowercase().contains("limits of usage") {
        return trimmed.to_string();
    }
    match trimmed.chars().last() {
        Some('.' | '!' | '?') => format!("{trimmed} Includes limits of usage."),
        _ => format!("{trimmed}. Includes limits of usage."),
    }
}

fn insert_data_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
    figure_count: usize,
    figure_manifest_path: &Path,
    asset_review_manifest_path: &Path,
    validation_rows_path: &Path,
) {
    let has_robustness = !input.option2_robustness.is_empty();
    let abstract_validation_sentence = if has_robustness {
        "Selected designs were screened for robustness under +/-10%/+/-20% operating-parameter perturbations.".to_string()
    } else {
        "No standalone perturbation-sweep robustness dataset was emitted for this canonical run."
            .to_string()
    };
    let methods_pipeline_step7 = if has_robustness {
        "7. **Robustness screening**: selected designs are ranked and screened under operating-parameter perturbations.".to_string()
    } else {
        String::new()
    };
    let milestone_scope_sentence = if has_robustness {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This report covers the CFDrs pipeline from candidate generation through strict gating, scoring, ranking, and robustness screening.".to_string()
    } else {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This report covers the CFDrs pipeline from candidate generation through strict gating, scoring, and ranking.".to_string()
    };
    values.insert(
        "TOTAL_CANDIDATES".to_string(),
        input.total_candidates.to_string(),
    );
    values.insert(
        "RUN_AUTHORITY_NOTE".to_string(),
        if input.fast_mode {
            format!(
                "Selected figures and ranked conclusions reflect the authoritative fast-mode subset of the Milestone 12 parameter lattice. The full canonical lattice spans {} candidates, but this report's ranked outputs are generated from deterministic stride-selected topology-family coverage plus dense fill inside the winning family rather than exhaustive evaluation of every lattice point.",
                input.total_candidates,
            )
        } else if input.authoritative_run {
            format!(
                "Selected figures and ranked conclusions reflect the complete Milestone 12 topology sweep across all {} candidates in the canonical parameter lattice.",
                input.total_candidates,
            )
        } else {
            format!(
                "This is a non-authoritative fast-mode run derived from a canonical parameter lattice of {} candidates. Option 1 topology coverage may be truncated relative to the full sweep.",
                input.total_candidates,
            )
        },
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
        abstract_validation_sentence,
    );
    values.insert("METHODS_PIPELINE_STEP7".to_string(), methods_pipeline_step7);
    values.insert(
        "MILESTONE_SCOPE_SENTENCE".to_string(),
        milestone_scope_sentence,
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
    values.insert(
        "VALIDATION_TRACEABILITY_INTRO".to_string(),
        build_validation_traceability_intro(
            &input.canonical_source,
            Some(figure_count),
            &figure_manifest_path.display().to_string(),
            &asset_review_manifest_path.display().to_string(),
            &validation_rows_path.display().to_string(),
            input.validation_rows,
        ),
    );
}

/// Wrap a markdown table in a centered div with a numbered caption.
fn centered_table(caption: &str, table: &str) -> String {
    let rendered_table = if looks_like_markdown_table(table) {
        markdown_table_to_compact_html(table)
    } else {
        table.to_string()
    };

    format!(
        "<div align=\"center\">\n\n<p><strong>{caption}</strong></p>\n\n{rendered_table}\n</div>"
    )
}

fn looks_like_markdown_table(table: &str) -> bool {
    let mut lines = table.lines().filter(|line| !line.trim().is_empty());
    matches!(lines.next(), Some(line) if line.trim_start().starts_with('|'))
        && matches!(lines.next(), Some(line) if line.contains("---"))
}

fn markdown_table_to_compact_html(table: &str) -> String {
    let mut lines = table.lines().filter(|line| !line.trim().is_empty());
    let Some(header_line) = lines.next() else {
        return table.to_string();
    };
    let Some(separator_line) = lines.next() else {
        return table.to_string();
    };
    if !separator_line.contains("---") {
        return table.to_string();
    }

    let headers = parse_markdown_table_row(header_line);
    if headers.is_empty() {
        return table.to_string();
    }

    let rows: Vec<Vec<String>> = lines
        .filter(|line| line.trim_start().starts_with('|'))
        .map(parse_markdown_table_row)
        .collect();

    let mut html = String::new();
    html.push_str(
        "<table style=\"width:88%; max-width:6.2in; margin:0 auto; border-collapse:collapse; table-layout:fixed; font-size:8.5pt; line-height:1.15;\">\n<thead><tr>",
    );
    for header in &headers {
        html.push_str(&format!(
            "<th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:{};\">{}</th>",
            column_width_for_header(header),
            render_table_cell_inline_markdown(header)
        ));
    }
    html.push_str("</tr></thead>\n<tbody>");

    for row in rows {
        html.push_str("<tr>");
        for (index, cell) in row.iter().enumerate() {
            let header = headers.get(index).map_or("", String::as_str);
            let alignment = if is_numeric_column(header) {
                "center"
            } else {
                "left"
            };
            html.push_str(&format!(
                "<td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:{}; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{}</td>",
                alignment,
                render_table_cell_inline_markdown(cell)
            ));
        }
        html.push_str("</tr>\n");
    }
    html.push_str("</tbody></table>");
    html
}

fn parse_markdown_table_row(row: &str) -> Vec<String> {
    row.trim()
        .trim_matches('|')
        .split('|')
        .map(|cell| cell.trim().to_string())
        .collect()
}

fn column_width_for_header(header: &str) -> &'static str {
    let normalized = header.to_ascii_lowercase();
    if normalized.contains("candidate") {
        "22%"
    } else if normalized.contains("track") || normalized.contains("topology") {
        "14%"
    } else if normalized.contains("mode") {
        "10%"
    } else if normalized.contains("parameter") {
        "18%"
    } else if normalized.contains("reference") || normalized.contains("model") {
        "16%"
    } else if normalized.contains("crate") {
        "10%"
    } else {
        "8%"
    }
}

fn is_numeric_column(header: &str) -> bool {
    let normalized = header.to_ascii_lowercase();
    normalized.contains("score")
        || normalized.contains("sigma")
        || normalized.contains("dose")
        || normalized.contains("exposure")
        || normalized.contains("fraction")
        || normalized.contains("efficiency")
        || normalized.contains("count")
        || normalized.contains("rank")
        || normalized.contains("shear")
        || normalized.contains("ecv")
        || normalized.contains("risk")
        || normalized.contains("cv")
        || normalized.contains("min")
        || normalized.contains("max")
        || normalized.contains("status")
}

fn render_table_cell_inline_markdown(text: &str) -> String {
    let bytes = text.as_bytes();
    let mut out = String::new();
    let mut idx = 0;
    while idx < bytes.len() {
        if idx + 1 < bytes.len() && bytes[idx] == b'*' && bytes[idx + 1] == b'*' {
            if let Some(end) = text[idx + 2..].find("**") {
                out.push_str("<strong>");
                out.push_str(&html_escape(&text[idx + 2..idx + 2 + end]));
                out.push_str("</strong>");
                idx += 2 + end + 2;
                continue;
            }
        }
        if bytes[idx] == b'`' {
            if let Some(end) = text[idx + 1..].find('`') {
                out.push_str("<code>");
                out.push_str(&html_escape(&text[idx + 1..idx + 1 + end]));
                out.push_str("</code>");
                idx += 1 + end + 1;
                continue;
            }
        }
        let ch = text[idx..].chars().next().expect("valid char boundary");
        match ch {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(ch),
        }
        idx += ch.len_utf8();
    }
    out
}

fn html_escape(text: &str) -> String {
    text.chars()
        .map(|ch| match ch {
            '&' => "&amp;".to_string(),
            '<' => "&lt;".to_string(),
            '>' => "&gt;".to_string(),
            '"' => "&quot;".to_string(),
            '\'' => "&#39;".to_string(),
            _ => ch.to_string(),
        })
        .collect()
}

fn insert_table_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) {
    let has_robustness = !input.option2_robustness.is_empty();
    let validation_table_number = if has_robustness { 14 } else { 13 };
    values.insert(
        "RESULTS_INTRO".to_string(),
        build_results_intro(
            input.total_candidates,
            input.option1_evaluated_count,
            input.option2_evaluated_count,
            input.option1_pool_len,
            input.option2_pool_len,
            input.fast_mode,
            has_robustness,
        ),
    );
    values.insert(
        "VALIDATION_SUMMARY_TABLE".to_string(),
        centered_table(
            &format!("Table {validation_table_number}. Cross-fidelity validation evidence"),
            &build_validation_summary_table(input.validation_rows),
        ),
    );
    values.insert(
        "SELECTED_TABLE_INTRO".to_string(),
        build_selected_table_intro(option1, option2),
    );
    values.insert(
        "SELECTED_TABLE".to_string(),
        centered_table(
            "Table 8. Selected designs comparison",
            &build_selected_table(option1, option2),
        ),
    );
    values.insert(
        "TOP5_INTRO".to_string(),
        build_top5_intro(input.option2_ranked),
    );
    values.insert(
        "OPTION2_TOP5_TABLE".to_string(),
        centered_table(
            "Table 9. Option 2 top-5 deterministic ranking",
            &build_option2_top5_table(input.option2_ranked),
        ),
    );
    values.insert(
        "OPTION1_SEQUENCE_FAMILY_SUMMARY".to_string(),
        centered_table(
            "Table 10. Option 1 topology-family coverage and ranking",
            &input.option1_sequence_summary_markdown,
        ),
    );
    values.insert(
        "TRI_CELL_INTRO".to_string(),
        build_tri_cell_intro(option2, input.ga_top.first()),
    );
    values.insert(
        "TRI_CELL_TABLE".to_string(),
        centered_table(
            "Table 11. Three-population routing evidence",
            &build_tri_cell_table(input.option2_ranked, input.ga_top),
        ),
    );
    values.insert("STRICT_CORE_INTRO".to_string(), build_strict_core_intro());
    values.insert(
        "STRICT_CORE_TABLE".to_string(),
        centered_table(
            "Table 12. Eligibility gate evidence",
            &build_strict_core_table(input.option2_ranked, input.ga_top),
        ),
    );
    values.insert(
        "ROBUSTNESS_SECTION".to_string(),
        build_robustness_section(input.option2_robustness, input.fast_mode, 13),
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
        "TREATMENT_TIME_ANALYSIS".to_string(),
        build_treatment_time_analysis(option2, if has_robustness { 16 } else { 15 }),
    );
    values.insert(
        "CRI_EXPANSION_SENSITIVITY".to_string(),
        build_cri_expansion_sensitivity(if has_robustness { 17 } else { 16 }),
    );
    values.insert(
        "PHYSICS_MODEL_INVENTORY".to_string(),
        physics_model_inventory(),
    );
    values.insert(
        "SERPENTINE_VENTURI_DISCUSSION".to_string(),
        build_serpentine_venturi_discussion(option2, input.ga_top.first()),
    );
}
