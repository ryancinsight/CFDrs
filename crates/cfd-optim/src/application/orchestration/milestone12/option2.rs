use std::collections::HashMap;
use std::fs;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::analysis::{robustness_sweep_blueprint, RobustnessReport, STANDARD_PERTURBATIONS};
use crate::application::objectives::{
    score_selective_acoustic_residence_separation, score_selective_venturi_cavitation,
};
use crate::application::orchestration::{
    blueprint_lineage_key, ensure_release_reports, fast_env, fast_mode, init_tracing,
    resolve_output_directories, save_figure, ScanProgress,
};
use crate::delivery::{load_top5_report_json, save_pareto_points, save_top5_report_json};
use crate::design::{build_milestone12_candidate_params, CandidateParams};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::metrics::evaluate_blueprint_candidate;
use crate::reporting::{
    audit_goal_candidates, validate_milestone12_candidate, write_goal_audit_report,
    Milestone12LineageKey, Milestone12ReportDesign, Milestone12Stage, ParetoPoint, ParetoTag,
};

use super::report::{write_stage_summary, Milestone12Option2Summary, OPTION2_SUMMARY_PATH};
use super::types::{Milestone12Option2Run, Milestone12StageArtifact};

/// Lightweight result stored during parallel evaluation — holds the candidate
/// ID and lineage key instead of the full ~15 KB `BlueprintCandidate`.
/// The full candidate is only re-materialized for the final top-5.
struct LightweightResult {
    params: CandidateParams,
    lineage_key: Option<Milestone12LineageKey>,
    is_venturi: bool,
    option2_score: Option<f64>,
}

pub fn run_milestone12_option2() -> Result<Milestone12Option2Run, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, figures_dir) = resolve_output_directories()?;
    let is_fast = fast_mode();
    let goal = OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity;

    // Keep lightweight params (~100 bytes each) instead of full candidates (~15 KB).
    // At 500K entries: ~50 MB vs ~7.5 GB.
    let all_params = build_milestone12_candidate_params();
    let total_candidates = all_params.len();
    let option1_from_disk =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let have_option1 = !option1_from_disk.is_empty();

    // In fast mode, cap evaluation to avoid hanging on the full space.
    let eval_cap = if is_fast {
        fast_env("M12_FAST_VENTURI_EVAL_MAX", 1000)
    } else {
        total_candidates
    };

    // ── Focused/strided optimization (3 phases) ──
    // Phase 1: Coarse stride across ALL params → find best topology family
    // Phase 2: Medium stride within winning family → find best stride pair
    // Phase 3: Dense fill between the two best adjacent stride points
    let all_params_vec = all_params;
    let phase1_budget = eval_cap / 4;
    let phase2_budget = eval_cap / 4;
    let phase3_budget = eval_cap / 2;

    // Phase 1: coarse stride
    let stride1 = (all_params_vec.len() / phase1_budget.max(1)).max(1);
    let phase1_params: Vec<CandidateParams> = all_params_vec
        .iter()
        .step_by(stride1)
        .take(phase1_budget)
        .cloned()
        .collect();

    // Audit on phase1 subset.
    let audit_cap = phase1_params.len().min(100);
    let selective_for_audit: Vec<BlueprintCandidate> = phase1_params
        .iter()
        .take(audit_cap)
        .map(|p| p.materialize())
        .collect();
    let audit_entries = audit_goal_candidates(&selective_for_audit, goal);
    let audit = write_goal_audit_report(&out_dir, "option2_audit", &audit_entries)?;
    drop(selective_for_audit);
    drop(audit_entries);

    // Phase 1 scoring: find best topology family.
    let mut best_score = f64::NEG_INFINITY;
    let mut winning_family = String::new();
    for p in &phase1_params {
        let candidate = p.materialize();
        if let Ok(eval) = evaluate_blueprint_candidate(&candidate) {
            let is_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
            if is_venturi {
                if let Some(score) = score_selective_venturi_cavitation(&eval, true) {
                    if score > best_score {
                        best_score = score;
                        winning_family = p.seq_tag();
                    }
                }
            }
        }
    }

    // Phase 2: stride within winning family
    let family_params: Vec<CandidateParams> = all_params_vec
        .iter()
        .filter(|p| p.seq_tag() == winning_family)
        .cloned()
        .collect();
    let stride2 = (family_params.len() / phase2_budget.max(1)).max(1);
    let phase2_indices: Vec<usize> = (0..family_params.len())
        .step_by(stride2)
        .take(phase2_budget)
        .collect();

    // Score phase2 stride points, find best pair.
    let mut scored_indices: Vec<(usize, f64)> = Vec::with_capacity(phase2_indices.len());
    for &i in &phase2_indices {
        let candidate = family_params[i].materialize();
        if let Ok(eval) = evaluate_blueprint_candidate(&candidate) {
            if let Some(score) = score_selective_venturi_cavitation(&eval, true) {
                scored_indices.push((i, score));
            }
        }
    }
    scored_indices.sort_by(|a, b| b.1.total_cmp(&a.1));
    let mut best_indices: Vec<usize> = scored_indices.iter().take(2).map(|&(i, _)| i).collect();
    best_indices.sort();

    // Phase 3: dense fill between best stride pair
    let (fill_lo, fill_hi) = if best_indices.len() >= 2 {
        (best_indices[0].saturating_sub(stride2), best_indices[1] + stride2)
    } else if let Some(&idx) = best_indices.first() {
        (idx.saturating_sub(stride2 * 2), (idx + stride2 * 2).min(family_params.len()))
    } else {
        (0, family_params.len().min(phase3_budget))
    };
    let phase3_params: Vec<CandidateParams> = family_params
        [fill_lo..fill_hi.min(family_params.len())]
        .iter()
        .take(phase3_budget)
        .cloned()
        .collect();

    // Merge all phases, deduplicate by idx.
    let mut seen_idx = std::collections::HashSet::with_capacity(eval_cap);
    let selective_params: Vec<CandidateParams> = phase1_params
        .into_iter()
        .chain(phase2_indices.iter().map(|&i| family_params[i].clone()))
        .chain(phase3_params)
        .filter(|p| seen_idx.insert(p.idx))
        .collect();

    let progress = Arc::new(ScanProgress::new("option2 scan", selective_params.len()));
    let have_option1_for_par = have_option1;

    // --- Streaming evaluation with lazy materialization ---
    // Each candidate is materialized, evaluated, and the full object dropped
    // inside the parallel closure. Only lightweight data is retained.
    let pareto_acc = Mutex::new(Vec::with_capacity(selective_params.len() / 2));
    let deferred_acc = Mutex::new(Vec::<LightweightResult>::new());

    selective_params
        .into_par_iter()
        .for_each(|params| {
            let progress_ref = &progress;
            let candidate = params.materialize();
            let evaluation = match evaluate_blueprint_candidate(&candidate) {
                Ok(evaluation) => evaluation,
                Err(_) => {
                    progress_ref.record();
                    return;
                }
            };
            progress_ref.record();

            let is_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);
            let lineage_key = blueprint_lineage_key(&candidate);

            // Score-only functions: borrow &evaluation, zero allocation, no clone.
            let option1_score = if !is_venturi && !have_option1_for_par {
                score_selective_acoustic_residence_separation(&evaluation)
            } else {
                None
            };
            let option2_score = if is_venturi {
                let has_placements = candidate
                    .topology_spec()
                    .map_or(false, |t| !t.venturi_placements.is_empty());
                score_selective_venturi_cavitation(&evaluation, has_placements)
            } else {
                None
            };

            // Compute lightweight Pareto data inline from evaluation.
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

            // Drop the full candidate + evaluation here (end of closure).
            // Only push lightweight data to accumulators.
            if is_venturi {
                pareto_acc.lock().unwrap().push(ParetoPoint {
                    cancer_targeted_cavitation,
                    rbc_venturi_protection,
                    score: option2_score.unwrap_or(0.001),
                    tag: ParetoTag::Option2,
                });
            }

            if option1_score.is_some() || option2_score.is_some() {
                deferred_acc.lock().unwrap().push(LightweightResult {
                    params,
                    lineage_key,
                    is_venturi,
                    option2_score,
                });
            }
            // `candidate` and `evaluation` are dropped here — no accumulation.
        });
    progress.finish();

    let pareto_points = pareto_acc.into_inner().unwrap();
    let deferred = deferred_acc.into_inner().unwrap();

    // Persist lightweight Pareto points (~32 bytes each).
    let option2_pool_all_path = out_dir.join("option2_pool_all.json");
    save_pareto_points(&pareto_points, &option2_pool_all_path)?;
    drop(pareto_points);

    // --- Lineage matching on lightweight data ---
    let mut option1_lineage_keys: HashMap<Milestone12LineageKey, bool> = HashMap::new();
    if have_option1 {
        for design in &option1_from_disk {
            if let Some(key) = blueprint_lineage_key(&design.candidate) {
                option1_lineage_keys.insert(key, true);
            }
        }
    } else {
        for result in &deferred {
            if !result.is_venturi {
                if let Some(ref key) = result.lineage_key {
                    option1_lineage_keys.insert(key.clone(), true);
                }
            }
        }
    }

    // Split into option1 and option2 deferred (lightweight — only params + score).
    let option2_deferred: Vec<(usize, f64)> = deferred
        .iter()
        .enumerate()
        .filter_map(|(i, r)| {
            if r.is_venturi {
                r.option2_score.map(|s| (i, s))
            } else {
                None
            }
        })
        .collect();

    let mut option2_by_lineage: HashMap<Milestone12LineageKey, Vec<usize>> = HashMap::new();
    for &(deferred_idx, _) in &option2_deferred {
        if let Some(ref key) = deferred[deferred_idx].lineage_key {
            option2_by_lineage.entry(key.clone()).or_default().push(deferred_idx);
        }
    }

    let mut viable_lineages: Vec<(Milestone12LineageKey, Vec<usize>)> = option2_by_lineage
        .into_iter()
        .filter(|(key, _)| option1_lineage_keys.contains_key(key))
        .collect();
    drop(option1_lineage_keys);

    // Select top-5 indices into the `deferred` vec.
    let selected_deferred_indices: Vec<usize> = if viable_lineages.is_empty() {
        // No shared lineage — sort all venturi by score, take top 5.
        let mut sorted: Vec<(usize, f64)> = option2_deferred;
        sorted.sort_by(|a, b| b.1.total_cmp(&a.1));
        sorted.into_iter().take(5).map(|(i, _)| i).collect()
    } else {
        // Sort viable lineages by best Option 2 score — O(n) direct index
        // lookup instead of O(n²) find_map scan.
        viable_lineages.sort_by(|left, right| {
            let left_best = left
                .1
                .iter()
                .map(|&i| deferred[i].option2_score.unwrap_or(f64::NEG_INFINITY))
                .fold(f64::NEG_INFINITY, f64::max);
            let right_best = right
                .1
                .iter()
                .map(|&i| deferred[i].option2_score.unwrap_or(f64::NEG_INFINITY))
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
        // Sort by score descending via lookup into deferred.
        winning.sort_by(|&a, &b| {
            let sa = deferred[a].option2_score.unwrap_or(f64::NEG_INFINITY);
            let sb = deferred[b].option2_score.unwrap_or(f64::NEG_INFINITY);
            sb.total_cmp(&sa)
        });
        winning.truncate(5);
        winning
    };
    drop(viable_lineages);

    // --- Re-materialize ONLY the final top-5 candidates ---
    // 5 × 15 KB = 75 KB instead of 50K × 15 KB = 750 MB.
    let mut option2_pool: Vec<Milestone12ReportDesign> =
        Vec::with_capacity(selected_deferred_indices.len());
    for &idx in &selected_deferred_indices {
        let result = &deferred[idx];
        let score = result.option2_score.unwrap_or(0.001);
        let candidate = result.params.materialize();
        if let Ok(design) =
            Milestone12ReportDesign::from_blueprint_candidate(0, candidate, score)
        {
            option2_pool.push(design);
        }
    }
    drop(deferred); // free all lightweight results

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
