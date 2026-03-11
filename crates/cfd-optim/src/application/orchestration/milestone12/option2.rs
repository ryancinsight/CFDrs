use std::collections::HashMap;
use std::fs;
use std::sync::Arc;

use rayon::prelude::*;

use crate::analysis::{robustness_sweep_blueprint, RobustnessReport, STANDARD_PERTURBATIONS};
use crate::application::objectives::{
    evaluate_selective_acoustic_residence_separation, evaluate_selective_venturi_cavitation,
};
use crate::application::orchestration::{
    blueprint_lineage_key, ensure_release_reports, fast_mode, init_tracing,
    is_selective_report_topology, resolve_output_directories, save_figure, ScanProgress,
};
use crate::delivery::{load_top5_report_json, save_pareto_points, save_top5_report_json};
use crate::design::build_milestone12_blueprint_candidate_space;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::metrics::evaluate_blueprint_candidate;
use crate::reporting::{
    audit_goal_candidates, validate_milestone12_candidate, write_goal_audit_report,
    Milestone12LineageKey, Milestone12ReportDesign, Milestone12Stage, ParetoPoint, ParetoTag,
};

use super::report::{write_stage_summary, Milestone12Option2Summary, OPTION2_SUMMARY_PATH};
use super::types::{Milestone12Option2Run, Milestone12StageArtifact};

struct CandidateResult {
    candidate: BlueprintCandidate,
    is_venturi: bool,
    option1_score: Option<f64>,
    option2_score: Option<f64>,
    // Lightweight Pareto data extracted inline — avoids computing full SdtMetrics
    // for every candidate (~1.5 KB + a redundant 1-D solve per entry).
    cancer_targeted_cavitation: f64,
    rbc_venturi_protection: f64,
}

pub fn run_milestone12_option2() -> Result<Milestone12Option2Run, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, figures_dir) = resolve_output_directories()?;
    let is_fast = fast_mode();
    let goal = OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity;

    let all_candidates = build_milestone12_blueprint_candidate_space()?;
    let total_candidates = all_candidates.len();
    let option1_from_disk =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let have_option1 = !option1_from_disk.is_empty();

    let selective_candidates: Vec<BlueprintCandidate> = all_candidates
        .into_iter()
        .filter(is_selective_report_topology)
        .collect();

    let audit_entries = audit_goal_candidates(&selective_candidates, goal);
    let audit = write_goal_audit_report(&out_dir, "option2_audit", &audit_entries)?;

    let progress = Arc::new(ScanProgress::new("option2 scan", selective_candidates.len()));
    let progress_ref = Arc::clone(&progress);
    let have_option1_for_par = have_option1;

    let results: Vec<CandidateResult> = selective_candidates
        .into_par_iter()
        .filter_map(move |candidate| {
            let evaluation = match evaluate_blueprint_candidate(&candidate) {
                Ok(evaluation) => evaluation,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            progress_ref.record();

            let is_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
            let option1_score = if !is_venturi && !have_option1_for_par {
                let evaluation =
                    evaluate_selective_acoustic_residence_separation(&candidate, evaluation.clone());
                evaluation.score.filter(|_| evaluation.is_eligible())
            } else {
                None
            };
            let option2_score = if is_venturi {
                evaluate_selective_venturi_cavitation(&candidate, evaluation.clone())
                    .ok()
                    .filter(|evaluation| evaluation.is_eligible())
                    .and_then(|evaluation| evaluation.score)
            } else {
                None
            };

            // Compute lightweight Pareto data inline from evaluation —
            // avoids materialising full SdtMetrics (~1.5 KB + redundant 1-D solve).
            let cav_potential = evaluation
                .venturi
                .placements
                .iter()
                .map(|p| (1.0 - p.cavitation_number).clamp(0.0, 1.0))
                .fold(0.0_f64, f64::max);
            let constriction = evaluation.venturi.cavitation_selectivity_score.clamp(0.0, 1.0);
            let cavitation_intensity = cav_potential * (0.5 + 0.5 * constriction);
            let cancer_targeted_cavitation =
                evaluation.separation.cancer_center_fraction.clamp(0.0, 1.0) * cavitation_intensity;
            let rbc_venturi_protection =
                evaluation.separation.rbc_peripheral_fraction.clamp(0.0, 1.0)
                    * (1.0
                        - cavitation_intensity
                            * evaluation.venturi.rbc_exposure_fraction.clamp(0.0, 1.0));

            Some(CandidateResult {
                candidate,
                is_venturi,
                option1_score,
                option2_score,
                cancer_targeted_cavitation,
                rbc_venturi_protection,
            })
        })
        .collect();
    progress.finish();

    // Separate results into lightweight deferred pools.
    // Full SdtMetrics are computed ONLY for the final top-5 (not all ~50K).
    let mut pareto_points = Vec::new();
    let mut option1_deferred: Vec<(BlueprintCandidate, f64)> = Vec::new();
    // (candidate, score) — no metrics yet
    let mut option2_deferred: Vec<(BlueprintCandidate, f64)> = Vec::new();
    for result in results {
        if result.is_venturi {
            pareto_points.push(ParetoPoint {
                cancer_targeted_cavitation: result.cancer_targeted_cavitation,
                rbc_venturi_protection: result.rbc_venturi_protection,
                score: result.option2_score.unwrap_or(0.001),
                tag: ParetoTag::Option2,
            });
            if let Some(option2_score) = result.option2_score {
                option2_deferred.push((result.candidate, option2_score));
            }
        } else if let Some(option1_score) = result.option1_score {
            option1_deferred.push((result.candidate, option1_score));
        }
    }

    // Persist lightweight Pareto points (~32 bytes each vs ~20 KB for full designs).
    let option2_pool_all_path = out_dir.join("option2_pool_all.json");
    save_pareto_points(&pareto_points, &option2_pool_all_path)?;
    drop(pareto_points);

    // --- Lineage matching on lightweight (candidate, score) pairs ---
    // No SdtMetrics needed — only topology spec keys and scores.
    let mut option1_lineage_keys: HashMap<Milestone12LineageKey, bool> = HashMap::new();
    if have_option1 {
        for design in &option1_from_disk {
            if let Some(key) = blueprint_lineage_key(&design.candidate) {
                option1_lineage_keys.insert(key, true);
            }
        }
    } else {
        for (candidate, _score) in &option1_deferred {
            if let Some(key) = blueprint_lineage_key(candidate) {
                option1_lineage_keys.insert(key, true);
            }
        }
    }
    drop(option1_deferred); // free option1 candidates early

    let mut option2_by_lineage: HashMap<Milestone12LineageKey, Vec<usize>> = HashMap::new();
    for (index, (candidate, _score)) in option2_deferred.iter().enumerate() {
        if let Some(key) = blueprint_lineage_key(candidate) {
            option2_by_lineage.entry(key).or_default().push(index);
        }
    }

    let mut viable_lineages: Vec<(Milestone12LineageKey, Vec<usize>)> = option2_by_lineage
        .into_iter()
        .filter(|(key, _)| option1_lineage_keys.contains_key(key))
        .collect();
    drop(option1_lineage_keys);

    // Select and truncate the lightweight pool — no metrics computed yet.
    let selected_indices: Vec<usize> = if viable_lineages.is_empty() {
        // No shared lineage — sort all by score, take top 5.
        option2_deferred.sort_by(|a, b| b.1.total_cmp(&a.1));
        (0..option2_deferred.len().min(5)).collect()
    } else {
        viable_lineages.sort_by(|left, right| {
            let left_best = left
                .1
                .iter()
                .map(|&i| option2_deferred[i].1)
                .fold(f64::NEG_INFINITY, f64::max);
            let right_best = right
                .1
                .iter()
                .map(|&i| option2_deferred[i].1)
                .fold(f64::NEG_INFINITY, f64::max);
            right_best.total_cmp(&left_best)
        });
        let mut winning = viable_lineages
            .iter()
            .find(|(_, indices)| indices.len() >= 5)
            .or(viable_lineages.first())
            .ok_or("no shared Option 1 → Option 2 lineage retained venturi candidates")?
            .1
            .clone();
        // Sort winning indices by score descending, take top 5.
        winning.sort_by(|&a, &b| option2_deferred[b].1.total_cmp(&option2_deferred[a].1));
        winning.truncate(5);
        winning
    };
    drop(viable_lineages);

    // Materialise full Milestone12ReportDesign ONLY for the final top-5.
    // This is where SdtMetrics (~1.5 KB + 1-D solve) are computed —
    // for 5 candidates instead of ~50K.
    let mut option2_pool: Vec<Milestone12ReportDesign> = Vec::with_capacity(selected_indices.len());
    // Move selected candidates out; swap_remove would invalidate indices,
    // so we drain by sorted index in reverse.
    let mut to_take: Vec<(usize, BlueprintCandidate, f64)> = {
        let mut idx_sorted = selected_indices.clone();
        idx_sorted.sort_unstable();
        idx_sorted.dedup();
        // Take from highest index first so earlier indices stay valid.
        let mut taken = Vec::with_capacity(idx_sorted.len());
        for &i in idx_sorted.iter().rev() {
            let (candidate, score) = option2_deferred.swap_remove(i);
            taken.push((i, candidate, score));
        }
        taken
    };
    // Restore original selection order.
    to_take.sort_by_key(|(original_idx, _, _)| {
        selected_indices.iter().position(|&si| si == *original_idx).unwrap_or(usize::MAX)
    });
    drop(option2_deferred); // free the bulk of deferred candidates

    for (_, candidate, score) in to_take {
        if let Ok(design) = Milestone12ReportDesign::from_blueprint_candidate(0, candidate, score) {
            option2_pool.push(design);
        }
    }

    for (index, design) in option2_pool.iter_mut().enumerate() {
        design.rank = index + 1;
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option2Derived)?;
    }
    let eligible_count = option2_pool.len();
    if option2_pool.is_empty() {
        return Err("No eligible Option 2 candidates".into());
    }

    let top5_path = out_dir.join("two_concept_option2_venturi_top5.json");
    save_top5_report_json(&option2_pool, &top5_path)?;

    let robustness_path = out_dir.join("option2_combined_robustness_top5.json");
    let robustness = if is_fast {
        Vec::new()
    } else {
        let robustness: Vec<RobustnessReport> = option2_pool
            .iter()
            .map(|design| {
                robustness_sweep_blueprint(
                    &design.candidate,
                    OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
                    &STANDARD_PERTURBATIONS,
                )
            })
            .collect();
        fs::write(&robustness_path, serde_json::to_string_pretty(&robustness)?)?;
        robustness
    };

    let figure_path = figures_dir.join("selected_option2_combined_schematic.svg");
    save_figure(
        option2_pool[0].candidate.blueprint(),
        &figure_path,
        "Figure 5 (Option 2 combined selective venturi)",
    );

    write_stage_summary(
        &out_dir,
        OPTION2_SUMMARY_PATH,
        &Milestone12Option2Summary {
            total_candidates,
            eligible_count,
            authoritative_run: !is_fast,
            fast_mode: is_fast,
        },
    )?;

    Ok(Milestone12Option2Run {
        total_candidates,
        eligible_count,
        ranked: option2_pool,
        robustness,
        audit,
        artifacts: vec![
            Milestone12StageArtifact {
                label: "option2_top5".to_string(),
                path: top5_path,
            },
            Milestone12StageArtifact {
                label: "option2_robustness".to_string(),
                path: robustness_path,
            },
            Milestone12StageArtifact {
                label: "option2_figure".to_string(),
                path: figure_path,
            },
        ],
    })
}
