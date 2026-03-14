use serde::Serialize;
use std::collections::{BTreeMap, HashMap};

use crate::application::orchestration::{
    ensure_release_reports, fast_env, fast_mode, init_tracing, is_selective_report_topology,
    resolve_output_directories, save_figure,
};
use crate::application::search::pool::EvaluatedPool;
use crate::delivery::{save_json_pretty, save_top5_report_json};
use crate::design::{build_milestone12_candidate_params, CandidateParams};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::reporting::{
    audit_goal_candidates, validate_milestone12_candidate, write_goal_audit_report,
    GoalAuditEntry, GoalAuditStatus, Milestone12ReportDesign, Milestone12Stage,
};

use super::report::{
    write_stage_summary, Milestone12Option1Summary, Milestone12SequenceCoverage,
    OPTION1_SUMMARY_PATH,
};
use super::types::{Milestone12Option1Run, Milestone12StageArtifact};

#[derive(Debug, Serialize)]
struct ConceptTopline {
    concept: String,
    mode: String,
    available: bool,
    status: String,
    candidate_id: String,
    score: Option<f64>,
    topology: String,
    treatment_zone_mode: String,
    active_venturi_throat_count: Option<usize>,
    cavitation_number: Option<f64>,
    cancer_center_fraction: Option<f64>,
    therapeutic_window_score: Option<f64>,
    hemolysis_index_per_pass: Option<f64>,
    wall_shear_p95_pa: Option<f64>,
    mean_residence_time_s: Option<f64>,
    treatment_zone_dwell_time_s: Option<f64>,
    total_ecv_ml: Option<f64>,
    therapy_channel_fraction: Option<f64>,
}

fn concept_topline(concept: &str, mode: &str, d: &Milestone12ReportDesign) -> ConceptTopline {
    ConceptTopline {
        concept: concept.to_string(),
        mode: mode.to_string(),
        available: true,
        status: "selected".to_string(),
        candidate_id: d.candidate.id.clone(),
        score: Some(d.score),
        topology: d.topology_short_code(),
        treatment_zone_mode: d.metrics.treatment_zone_mode.clone(),
        active_venturi_throat_count: Some(d.metrics.active_venturi_throat_count),
        cavitation_number: Some(d.metrics.cavitation_number),
        cancer_center_fraction: Some(d.metrics.cancer_center_fraction),
        therapeutic_window_score: Some(d.metrics.therapeutic_window_score),
        hemolysis_index_per_pass: Some(d.metrics.hemolysis_index_per_pass),
        wall_shear_p95_pa: Some(d.metrics.wall_shear_p95_pa),
        mean_residence_time_s: Some(d.metrics.mean_residence_time_s),
        treatment_zone_dwell_time_s: Some(d.metrics.treatment_zone_dwell_time_s),
        total_ecv_ml: Some(d.metrics.total_ecv_ml),
        therapy_channel_fraction: Some(d.metrics.therapy_channel_fraction),
    }
}

#[derive(Default)]
struct SequenceCoverageAccumulator {
    total_candidates: usize,
    eligible_count: usize,
    top_candidate_id: Option<String>,
    top_score: Option<f64>,
    best_treatment_residence_time_s: Option<f64>,
    best_separation_efficiency: Option<f64>,
    limiter_counts: BTreeMap<String, usize>,
}

fn candidate_sequence_label(candidate: &BlueprintCandidate) -> String {
    candidate
        .topology_spec()
        .map_or_else(|_| "Unknown".to_string(), |spec| spec.stage_sequence_label())
}

fn summarize_sequence_coverage(
    candidates: &[BlueprintCandidate],
    audit_entries: &[GoalAuditEntry],
) -> Vec<Milestone12SequenceCoverage> {
    let mut coverage: BTreeMap<String, SequenceCoverageAccumulator> = BTreeMap::new();
    let mut candidate_sequences = BTreeMap::new();
    for candidate in candidates {
        let sequence_label = candidate_sequence_label(candidate);
        candidate_sequences.insert(candidate.id.clone(), sequence_label.clone());
        coverage.entry(sequence_label).or_default().total_candidates += 1;
    }

    for entry in audit_entries {
        let Some(sequence_label) = candidate_sequences.get(&entry.candidate_id) else {
            continue;
        };
        let accumulator = coverage.entry(sequence_label.clone()).or_default();
        match entry.status {
            GoalAuditStatus::Eligible => {
                accumulator.eligible_count += 1;
                if entry.score > accumulator.top_score {
                    accumulator.top_candidate_id = Some(entry.candidate_id.clone());
                    accumulator.top_score = entry.score;
                    accumulator.best_treatment_residence_time_s = entry.treatment_residence_time_s;
                    accumulator.best_separation_efficiency = entry.separation_efficiency;
                }
            }
            GoalAuditStatus::ScreenedOut | GoalAuditStatus::Errored => {
                for reason in &entry.reasons {
                    *accumulator.limiter_counts.entry(reason.clone()).or_insert(0) += 1;
                }
            }
        }
    }

    coverage
        .into_iter()
        .map(|(sequence_label, accumulator)| {
            let dominant_limiter = accumulator
                .limiter_counts
                .into_iter()
                .max_by_key(|(_, count)| *count)
                .map(|(reason, _)| reason)
                .unwrap_or_else(|| "eligible under current physics".to_string());
            Milestone12SequenceCoverage {
                sequence_label,
                total_candidates: accumulator.total_candidates,
                eligible_count: accumulator.eligible_count,
                top_candidate_id: accumulator.top_candidate_id,
                top_score: accumulator.top_score,
                best_treatment_residence_time_s: accumulator.best_treatment_residence_time_s,
                best_separation_efficiency: accumulator.best_separation_efficiency,
                dominant_limiter,
            }
        })
        .collect()
}

pub fn run_milestone12_option1() -> Result<Milestone12Option1Run, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, figures_dir) = resolve_output_directories()?;
    let is_fast = fast_mode();
    let goal = OptimizationGoal::AsymmetricSplitResidenceSeparation;

    // Keep lightweight params (~100 bytes each) instead of full candidates (~15 KB).
    let all_params = build_milestone12_candidate_params();
    let total_candidates = all_params.len();

    let eval_cap = if is_fast {
        fast_env("M12_FAST_ACOUSTIC_EVAL_MAX", 200)
    } else {
        total_candidates
    };

    // Filter to non-venturi selective topologies. Materialize transiently
    // for the filter check, then keep only the lightweight params.
    let selective_params: Vec<CandidateParams> = all_params
        .into_iter()
        .filter(|params| {
            // Quick venturi check from params — avoids full materialization.
            if params.is_venturi() {
                return false;
            }
            let candidate = params.materialize();
            is_selective_report_topology(&candidate)
        })
        .take(eval_cap)
        .collect();

    // Audit + sequence coverage need full candidates, but only transiently.
    let to_evaluate: Vec<BlueprintCandidate> =
        selective_params.iter().map(|p| p.materialize()).collect();
    let audit_entries = audit_goal_candidates(&to_evaluate, goal);
    let sequence_coverage = summarize_sequence_coverage(&to_evaluate, &audit_entries);
    let audit = write_goal_audit_report(&out_dir, "option1_audit", &audit_entries)?;

    let pool = EvaluatedPool::from_candidates(&to_evaluate);
    let shortlist_size = if is_fast { 3 } else { 5 };
    let eligible_count = pool.count_eligible(goal);
    let top_results = pool.top_k(shortlist_size, goal)?;
    drop(pool);

    // Consume `to_evaluate` — extract only the ≤5 candidates needed for the shortlist.
    let needed_ids: Vec<String> = top_results.iter().map(|e| e.candidate_id.clone()).collect();
    let mut candidate_map: HashMap<String, BlueprintCandidate> =
        HashMap::with_capacity(needed_ids.len());
    for candidate in to_evaluate {
        if needed_ids.iter().any(|id| *id == candidate.id) {
            candidate_map.insert(candidate.id.clone(), candidate);
        }
    }

    let mut ranked: Vec<Milestone12ReportDesign> = top_results
        .into_iter()
        .enumerate()
        .filter_map(|(index, evaluation)| {
            let candidate = candidate_map.remove(&evaluation.candidate_id)?;
            Milestone12ReportDesign::from_blueprint_candidate(
                index + 1,
                candidate,
                evaluation.score?,
            )
                .ok()
        })
        .collect();
    for (index, design) in ranked.iter_mut().enumerate() {
        design.rank = index + 1;
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option1Base)?;
    }

    let top5_path = out_dir.join("two_concept_option1_ultrasound_top5.json");
    save_top5_report_json(&ranked, &top5_path)?;
    let topline_path = out_dir.join("option1_concept_topline.json");
    if let Some(topline) = ranked.first().map(|design| {
        concept_topline(
            "Option 1: Selective acoustic center treatment",
            "AsymmetricSplitResidenceSeparation",
            design,
        )
    }) {
        save_json_pretty(&vec![topline], &topline_path)?;
    }

    let figure_path = figures_dir.join("selected_option1_schematic.svg");
    if let Some(best) = ranked.first() {
        save_figure(
            best.candidate.blueprint(),
            &figure_path,
            "Figure 4 (Option 1 selective acoustic)",
        );
    }

    write_stage_summary(
        &out_dir,
        OPTION1_SUMMARY_PATH,
        &Milestone12Option1Summary {
            total_candidates,
            eligible_count,
            authoritative_run: !is_fast,
            fast_mode: is_fast,
            sequence_coverage,
        },
    )?;

    let artifacts = vec![
        Milestone12StageArtifact {
            label: "option1_top5".to_string(),
            path: top5_path,
        },
        Milestone12StageArtifact {
            label: "option1_topline".to_string(),
            path: topline_path,
        },
        Milestone12StageArtifact {
            label: "option1_figure".to_string(),
            path: figure_path,
        },
    ];

    Ok(Milestone12Option1Run {
        total_candidates,
        eligible_count,
        ranked,
        audit,
        artifacts,
    })
}
