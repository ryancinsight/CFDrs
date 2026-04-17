use std::fs;
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::application::orchestration::{
    ensure_release_reports, fast_mode, init_tracing, resolve_output_directories,
};
use crate::delivery::{load_pareto_points, load_top5_report_json};
use crate::reporting::{
    pareto_pool_from_report_designs, rank_ga_hydrosdt_report_designs,
    write_milestone12_narrative_report, write_milestone12_results, Milestone12GaRankingAuditEntry,
    Milestone12NarrativeInput, ParetoTag, ValidationRow,
};

use super::ga::run_milestone12_ga;
use super::option1::run_milestone12_option1;
use super::option2::run_milestone12_option2;
use super::types::{Milestone12RequestedStage, Milestone12StageArtifact};
use super::validation::run_milestone12_validation;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub(crate) struct Milestone12SequenceCoverage {
    pub sequence_label: String,
    pub total_candidates: usize,
    pub eligible_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_candidate_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_score: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_treatment_residence_time_s: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_separation_efficiency: Option<f64>,
    #[serde(default)]
    pub dominant_limiter: String,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub(crate) struct Milestone12Option1Summary {
    pub total_candidates: usize,
    #[serde(default)]
    pub evaluated_count: usize,
    pub eligible_count: usize,
    #[serde(default)]
    pub authoritative_run: bool,
    #[serde(default)]
    pub fast_mode: bool,
    #[serde(default)]
    pub sequence_coverage: Vec<Milestone12SequenceCoverage>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub(crate) struct Milestone12Option2Summary {
    pub total_candidates: usize,
    #[serde(default)]
    pub evaluated_count: usize,
    pub eligible_count: usize,
    #[serde(default)]
    pub authoritative_run: bool,
    #[serde(default)]
    pub fast_mode: bool,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub(crate) struct Milestone12GaSummary {
    #[serde(default)]
    pub evaluated_count: usize,
    #[serde(default)]
    pub authoritative_run: bool,
    #[serde(default)]
    pub fast_mode: bool,
    #[serde(default)]
    pub best_per_generation: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct Milestone12ReportManifest {
    #[serde(default)]
    authoritative_run: bool,
    #[serde(default)]
    run_class: String,
    #[serde(default)]
    fast_mode: bool,
    #[serde(default)]
    review_complete: bool,
    #[serde(default)]
    canonical_source: String,
    #[serde(default)]
    requested_stages: Vec<String>,
    #[serde(default)]
    artifacts: Vec<Milestone12StageArtifact>,
}

pub(crate) const OPTION1_SUMMARY_PATH: &str = "option1_stage_summary.json";
pub(crate) const OPTION2_SUMMARY_PATH: &str = "option2_stage_summary.json";
pub(crate) const GA_SUMMARY_PATH: &str = "ga_stage_summary.json";
const REPORT_MANIFEST_PATH: &str = "report_manifest.json";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Milestone12RunClass {
    FastRefresh,
    AuthoritativeFastSubset,
    AuthoritativeFull,
}

impl Milestone12RunClass {
    const fn is_authoritative(self) -> bool {
        !matches!(self, Self::FastRefresh)
    }

    const fn requires_review(self) -> bool {
        self.is_authoritative()
    }

    const fn label(self) -> &'static str {
        match self {
            Self::FastRefresh => "fast_refresh",
            Self::AuthoritativeFastSubset => "authoritative_fast_subset",
            Self::AuthoritativeFull => "authoritative_full",
        }
    }
}

fn classify_run(
    requested_stages: &[Milestone12RequestedStage],
    is_fast_mode: bool,
) -> Milestone12RunClass {
    let refresh_only = requested_stages
        .iter()
        .all(|stage| matches!(stage, Milestone12RequestedStage::Refresh));
    if refresh_only {
        Milestone12RunClass::FastRefresh
    } else if is_fast_mode {
        Milestone12RunClass::AuthoritativeFastSubset
    } else {
        Milestone12RunClass::AuthoritativeFull
    }
}

fn canonical_results_path(workspace_root: &Path) -> std::path::PathBuf {
    workspace_root.join("report").join("milestone12_results.md")
}

fn canonical_source_label() -> String {
    "report/milestone12_results.md".to_string()
}

fn render_option1_sequence_coverage(summary: &Milestone12Option1Summary) -> String {
    let mut markdown = String::new();
    if summary.sequence_coverage.is_empty() {
        markdown.push_str("_No Option 1 sequence-coverage ledger was generated for this run._\n");
        return markdown;
    }

    markdown.push_str(
        "| Sequence | Candidates scanned | Eligible | Top candidate | Top score | Best residence (s) | Best separation | Dominant limiter |\n",
    );
    markdown.push_str("|---|---:|---:|---|---:|---:|---:|---|\n");
    for entry in &summary.sequence_coverage {
        markdown.push_str(&format!(
            "| {} | {} | {} | {} | {} | {} | {} | {} |\n",
            entry.sequence_label,
            entry.total_candidates,
            entry.eligible_count,
            entry
                .top_candidate_id
                .as_deref()
                .map_or("—".to_string(), |id| format!("`{id}`")),
            entry
                .top_score
                .map_or("—".to_string(), |score| format!("{score:.4}")),
            entry
                .best_treatment_residence_time_s
                .map_or("—".to_string(), |value| format!("{value:.4}")),
            entry
                .best_separation_efficiency
                .map_or("—".to_string(), |value| format!("{value:.4}")),
            if entry.dominant_limiter.is_empty() {
                "—"
            } else {
                entry.dominant_limiter.as_str()
            },
        ));
    }
    markdown
}

pub(crate) fn write_stage_summary<T: Serialize>(
    out_dir: &Path,
    name: &str,
    summary: &T,
) -> Result<(), Box<dyn std::error::Error>> {
    fs::write(out_dir.join(name), serde_json::to_string_pretty(summary)?)?;
    Ok(())
}

fn read_summary<T>(path: &Path) -> Result<T, Box<dyn std::error::Error>>
where
    T: for<'de> Deserialize<'de> + Default,
{
    if !path.exists() {
        return Ok(T::default());
    }
    Ok(serde_json::from_str(&fs::read_to_string(path)?)?)
}

fn read_validation_rows(out_dir: &Path) -> Result<Vec<ValidationRow>, Box<dyn std::error::Error>> {
    let path = out_dir.join("milestone12_validation_rows.json");
    if !path.exists() {
        return Ok(Vec::new());
    }
    Ok(serde_json::from_str(&fs::read_to_string(path)?)?)
}

fn read_option2_robustness(
    out_dir: &Path,
) -> Result<Vec<crate::analysis::RobustnessReport>, Box<dyn std::error::Error>> {
    let path = out_dir.join("option2_combined_robustness_top5.json");
    if !path.exists() {
        return Ok(Vec::new());
    }
    Ok(serde_json::from_str(&fs::read_to_string(path)?)?)
}

fn read_ga_ranking_audit(
    out_dir: &Path,
) -> Result<Vec<Milestone12GaRankingAuditEntry>, Box<dyn std::error::Error>> {
    let path = out_dir.join("ga_lineage_audit_top5.json");
    if !path.exists() {
        return Ok(Vec::new());
    }
    Ok(serde_json::from_str(&fs::read_to_string(path)?)?)
}

pub fn refresh_milestone12_reports(
    requested_stages: &[Milestone12RequestedStage],
) -> Result<Vec<Milestone12StageArtifact>, Box<dyn std::error::Error>> {
    let (workspace_root, out_dir, _) = resolve_output_directories()?;
    let is_fast_mode = fast_mode();
    let run_class = classify_run(requested_stages, is_fast_mode);

    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let ga_ranked = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;
    if option1_ranked.is_empty() {
        return Err("Milestone 12 report requires a non-empty Option 1 ranked pool".into());
    }
    if option2_ranked.is_empty() {
        return Err("Milestone 12 report requires a non-empty Option 2 ranked pool".into());
    }
    if ga_ranked.is_empty() {
        return Err("Milestone 12 report requires a non-empty GA ranked pool".into());
    }

    // Load lightweight Pareto points (~32 bytes each) instead of full
    // Milestone12ReportDesign (~20 KB each). At 50K candidates this saves ~1 GB.
    let option2_pool_all = load_pareto_points(&out_dir.join("option2_pool_all.json"))?;
    let ga_pool_all = load_pareto_points(&out_dir.join("ga_pool_all.json"))?;
    let ga_ranked_for_report = rank_ga_hydrosdt_report_designs(&ga_ranked, &option2_ranked[0]);
    let selected_ga_ranked = if ga_ranked_for_report.is_empty() {
        ga_ranked.clone()
    } else {
        ga_ranked_for_report
    };
    let ga_pool_for_report = if ga_pool_all.is_empty() {
        pareto_pool_from_report_designs(&selected_ga_ranked, ParetoTag::Ga, 200)
    } else {
        ga_pool_all
    };

    let option1_summary: Milestone12Option1Summary =
        read_summary(&out_dir.join(OPTION1_SUMMARY_PATH))?;
    let option2_summary: Milestone12Option2Summary =
        read_summary(&out_dir.join(OPTION2_SUMMARY_PATH))?;
    let ga_summary: Milestone12GaSummary = read_summary(&out_dir.join(GA_SUMMARY_PATH))?;
    let validation_rows = read_validation_rows(&out_dir)?;
    let option2_robustness = read_option2_robustness(&out_dir)?;
    let ga_ranking_audit = read_ga_ranking_audit(&out_dir)?;

    let canonical = canonical_results_path(&workspace_root);
    if let Some(parent) = canonical.parent() {
        fs::create_dir_all(parent)?;
    }
    write_milestone12_results(
        option1_summary
            .total_candidates
            .max(option2_summary.total_candidates),
        option1_summary.eligible_count.max(option1_ranked.len()),
        option2_summary.eligible_count.max(option2_ranked.len()),
        &option1_ranked,
        &option2_ranked,
        &selected_ga_ranked[0],
        &validation_rows,
        &option2_robustness,
        false,
        &canonical_source_label(),
        &canonical,
    )?;

    let mut narrative = write_milestone12_narrative_report(
        &workspace_root,
        &canonical,
        &Milestone12NarrativeInput {
            authoritative_run: false,
            canonical_source: canonical_source_label(),
            total_candidates: option1_summary
                .total_candidates
                .max(option2_summary.total_candidates),
            option1_evaluated_count: option1_summary.evaluated_count,
            option2_evaluated_count: option2_summary.evaluated_count,
            option1_pool_len: option1_summary.eligible_count.max(option1_ranked.len()),
            option2_pool_len: option2_summary.eligible_count.max(option2_ranked.len()),
            option1_sequence_summary_markdown: render_option1_sequence_coverage(&option1_summary),
            option1_ranked: &option1_ranked,
            option2_ranked: &option2_ranked,
            ga_top: &selected_ga_ranked,
            option2_pool_all: &option2_pool_all,
            ga_pool_all: &ga_pool_for_report,
            validation_rows: &validation_rows,
            option2_robustness: &option2_robustness,
            ga_best_per_gen: &ga_summary.best_per_generation,
            ga_ranking_audit: &ga_ranking_audit,
            topology_family_count: option1_summary.sequence_coverage.len().max(1),
            fast_mode: is_fast_mode,
        },
    )?;
    if run_class.is_authoritative() && narrative.asset_review_complete {
        narrative = write_milestone12_narrative_report(
            &workspace_root,
            &canonical,
            &Milestone12NarrativeInput {
                authoritative_run: true,
                canonical_source: canonical_source_label(),
                total_candidates: option1_summary
                    .total_candidates
                    .max(option2_summary.total_candidates),
                option1_evaluated_count: option1_summary.evaluated_count,
                option2_evaluated_count: option2_summary.evaluated_count,
                option1_pool_len: option1_summary.eligible_count.max(option1_ranked.len()),
                option2_pool_len: option2_summary.eligible_count.max(option2_ranked.len()),
                option1_sequence_summary_markdown: render_option1_sequence_coverage(
                    &option1_summary,
                ),
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_ranked,
                ga_top: &selected_ga_ranked,
                option2_pool_all: &option2_pool_all,
                ga_pool_all: &ga_pool_for_report,
                validation_rows: &validation_rows,
                option2_robustness: &option2_robustness,
                ga_best_per_gen: &ga_summary.best_per_generation,
                ga_ranking_audit: &ga_ranking_audit,
                topology_family_count: option1_summary.sequence_coverage.len().max(1),
                fast_mode: is_fast_mode,
            },
        )?;
    }
    let report_is_authoritative = run_class.is_authoritative() && narrative.asset_review_complete;
    if report_is_authoritative {
        write_milestone12_results(
            option1_summary
                .total_candidates
                .max(option2_summary.total_candidates),
            option1_summary.eligible_count.max(option1_ranked.len()),
            option2_summary.eligible_count.max(option2_ranked.len()),
            &option1_ranked,
            &option2_ranked,
            &selected_ga_ranked[0],
            &validation_rows,
            &option2_robustness,
            true,
            &canonical_source_label(),
            &canonical,
        )?;
    }

    let artifacts = vec![
        Milestone12StageArtifact {
            label: "canonical_results".to_string(),
            path: canonical,
        },
        Milestone12StageArtifact {
            label: "narrative_report".to_string(),
            path: narrative.narrative_path,
        },
        Milestone12StageArtifact {
            label: "figure_manifest".to_string(),
            path: narrative.figure_manifest_path,
        },
        Milestone12StageArtifact {
            label: "asset_review_manifest".to_string(),
            path: narrative.asset_review_manifest_path.clone(),
        },
    ];

    fs::write(
        out_dir.join(REPORT_MANIFEST_PATH),
        serde_json::to_string_pretty(&Milestone12ReportManifest {
            authoritative_run: report_is_authoritative,
            run_class: run_class.label().to_string(),
            fast_mode: is_fast_mode,
            review_complete: narrative.asset_review_complete,
            canonical_source: canonical_source_label(),
            requested_stages: requested_stages
                .iter()
                .map(|stage| stage.label().to_string())
                .collect(),
            artifacts: artifacts.clone(),
        })?,
    )?;

    if run_class.requires_review() && !narrative.asset_review_complete {
        return Err(format!(
            "Milestone 12 authoritative report review is incomplete; open and review each generated asset listed in {} before rerunning",
            narrative.asset_review_manifest_path.display()
        )
        .into());
    }

    Ok(artifacts)
}

pub fn run_milestone12_report(
    requested_stages: &[Milestone12RequestedStage],
) -> Result<Vec<Milestone12StageArtifact>, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    for stage in requested_stages {
        match stage {
            Milestone12RequestedStage::Option1 => {
                run_milestone12_option1()?;
            }
            Milestone12RequestedStage::Option2 => {
                run_milestone12_option2()?;
            }
            Milestone12RequestedStage::Ga => {
                run_milestone12_ga()?;
            }
            Milestone12RequestedStage::Validation => {
                run_milestone12_validation()?;
            }
            Milestone12RequestedStage::Refresh => {}
        }
    }

    refresh_milestone12_reports(requested_stages)
}

#[cfg(test)]
mod tests {
    use super::{classify_run, Milestone12RequestedStage, Milestone12RunClass};

    #[test]
    fn refresh_only_runs_are_non_authoritative() {
        let class = classify_run(&[Milestone12RequestedStage::Refresh], true);
        assert_eq!(class, Milestone12RunClass::FastRefresh);
        assert!(!class.is_authoritative());
    }

    #[test]
    fn non_refresh_fast_runs_are_authoritative_fast_subset() {
        let class = classify_run(&[Milestone12RequestedStage::Option1], true);
        assert_eq!(class, Milestone12RunClass::AuthoritativeFastSubset);
        assert!(class.is_authoritative());
        assert!(class.requires_review());
    }
}
