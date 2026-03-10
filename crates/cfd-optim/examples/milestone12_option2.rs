//! Milestone 12 — Option 2: Selective Venturi Cavitation.
//!
//! Evaluates venturi-equipped selective-routing topologies, ranks them by
//! the Milestone 12 venturi-cavitation selectivity score, and applies
//! robustness screening.
//!
//! The lineage constraint (shared topology family with Option 1) is enforced
//! when Option 1 results exist on disk; otherwise the best venturi family is
//! selected unconstrainedly.
//!
//! Produces:
//! - `report/milestone12/two_concept_option2_venturi_top5.json`
//! - `report/milestone12/option2_combined_robustness_top5.json`
//! - `report/milestone12/option2_scored_seeds.json`  (for GA warm-start)
//! - `report/figures/selected_option2_combined_schematic.svg`
//! - Regenerates the narrative report.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_option2 --no-default-features
//! M12_FAST=1 cargo run -p cfd-optim --example milestone12_option2 --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, compute_blueprint_report_metrics,
    evaluate_blueprint_candidate, evaluate_blueprint_genetic_refinement,
    evaluate_selective_acoustic_residence_separation, fast_mode, init_tracing,
    is_selective_report_topology, load_top5_report_json, option2_mode,
    orchestration_lineage_key as blueprint_lineage_key, report_eligible_venturi_oncology,
    resolve_output_directories, robustness_sweep_blueprint, save_figure, save_json_pretty,
    save_top5_report_json, score_candidate, sort_report_designs, validate_milestone12_candidate,
    write_milestone12_narrative_report, write_milestone12_results, BlueprintCandidate,
    Milestone12LineageKey, Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage,
    OptimizationGoal, RobustnessReport, ScanProgress, SdtMetrics, SdtWeights,
    STANDARD_PERTURBATIONS,
};
use rayon::prelude::*;
use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    init_tracing();
    let (workspace_root, out_dir, figures_dir) = resolve_output_directories()?;

    tracing::info!("=== Milestone 12 — Option 2: Selective Venturi Cavitation ===\n");
    let overall_start = Instant::now();
    let weights = SdtWeights::default();
    let oncology_mode = option2_mode();
    let is_fast = fast_mode();

    // ── Build candidate space ────────────────────────────────────────────────
    let t_build = Instant::now();
    let all_candidates = build_milestone12_blueprint_candidate_space()?;
    let total_candidates = all_candidates.len();
    tracing::info!(
        "Candidate space: {} blueprints ({:.1}s)",
        total_candidates,
        t_build.elapsed().as_secs_f64()
    );

    // ── Load Option 1 results (from prior run) for lineage matching ─────────
    let option1_from_disk =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let have_option1 = !option1_from_disk.is_empty();
    if have_option1 {
        tracing::info!(
            "Loaded {} Option 1 designs from disk for lineage matching",
            option1_from_disk.len()
        );
    } else {
        tracing::info!("No Option 1 results on disk — will evaluate acoustic candidates inline");
    }

    // ── Evaluate candidates ──────────────────────────────────────────────────
    // We need both venturi candidates (for Option 2) and acoustic candidates
    // (for lineage matching if Option 1 wasn't run separately).
    let selective_candidates: Vec<BlueprintCandidate> = all_candidates
        .into_iter()
        .filter(|c| is_selective_report_topology(c))
        .collect();

    let progress = Arc::new(ScanProgress::new(
        "option2 scan",
        selective_candidates.len(),
    ));
    let progress_ref = Arc::clone(&progress);
    let have_option1_for_par = have_option1;

    struct CandidateResult {
        candidate: BlueprintCandidate,
        metrics: SdtMetrics,
        is_venturi: bool,
        ga_score: f64,
        option1_score: Option<f64>,
        option2_score: Option<f64>,
    }

    let results: Vec<CandidateResult> = selective_candidates
        .into_par_iter()
        .filter_map(move |candidate| {
            let eval = match evaluate_blueprint_candidate(&candidate) {
                Ok(e) => e,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            let metrics = match compute_blueprint_report_metrics(&candidate) {
                Ok(m) => m,
                Err(_) => {
                    progress_ref.record();
                    return None;
                }
            };
            progress_ref.record();

            let is_venturi = candidate
                .topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi);

            // For non-venturi candidates that need option1 scoring, clone the
            // eval before consuming it in the GA path.
            let option1_score = if !is_venturi && !have_option1_for_par {
                let opt1 =
                    evaluate_selective_acoustic_residence_separation(&candidate, eval.clone());
                if opt1.score > 0.0 {
                    Some(opt1.score)
                } else {
                    None
                }
            } else {
                None
            };

            let ga_score =
                evaluate_blueprint_genetic_refinement(&candidate, eval).map_or(0.0, |e| e.score);

            let option2_score = if is_venturi && report_eligible_venturi_oncology(&metrics) {
                Some(score_candidate(&metrics, oncology_mode, &weights))
            } else {
                None
            };

            Some(CandidateResult {
                candidate,
                metrics,
                is_venturi,
                ga_score,
                option1_score,
                option2_score,
            })
        })
        .collect();
    progress.finish();

    // ── Split into pools ─────────────────────────────────────────────────────
    let mut option1_pool: Vec<Milestone12ReportDesign> = Vec::new();
    let mut option2_pool: Vec<Milestone12ReportDesign> = Vec::new();
    let mut scored_seeds: Vec<(f64, BlueprintCandidate)> = Vec::new();

    // Move candidates instead of cloning when only one destination needs
    // them.  Clone only when a venturi candidate qualifies for BOTH
    // scored_seeds AND option2_pool.
    for r in results {
        if r.is_venturi {
            let need_seed = r.ga_score > 0.0;
            match (need_seed, r.option2_score) {
                (true, Some(score)) => {
                    scored_seeds.push((r.ga_score, r.candidate.clone()));
                    option2_pool.push(Milestone12ReportDesign::new(
                        0,
                        r.candidate,
                        r.metrics,
                        score,
                    ));
                }
                (true, None) => {
                    scored_seeds.push((r.ga_score, r.candidate));
                }
                (false, Some(score)) => {
                    option2_pool.push(Milestone12ReportDesign::new(
                        0,
                        r.candidate,
                        r.metrics,
                        score,
                    ));
                }
                (false, None) => {}
            }
        } else if let Some(score) = r.option1_score {
            option1_pool.push(Milestone12ReportDesign::new(
                0,
                r.candidate,
                r.metrics,
                score,
            ));
        }
    }

    tracing::info!(
        "Option 2 eligible: {}  (acoustic inline: {})",
        option2_pool.len(),
        option1_pool.len()
    );

    // ── Lineage matching ─────────────────────────────────────────────────────
    // Use loaded Option 1 results if available, otherwise use inline-evaluated acoustic pool.
    let option1_for_lineage: &[Milestone12ReportDesign] = if have_option1 {
        &option1_from_disk
    } else {
        &option1_pool
    };

    // Index-based lineage matching — store indices into existing vecs
    // instead of cloning entire Milestone12ReportDesign objects (each ~10-20 KB).
    let mut option1_lineage_keys: HashMap<Milestone12LineageKey, bool> = HashMap::new();
    for design in option1_for_lineage {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option1_lineage_keys.insert(key, true);
        }
    }

    let mut option2_by_lineage: HashMap<Milestone12LineageKey, Vec<usize>> = HashMap::new();
    for (idx, design) in option2_pool.iter().enumerate() {
        if let Some(key) = blueprint_lineage_key(&design.candidate) {
            option2_by_lineage.entry(key).or_default().push(idx);
        }
    }

    // Find lineage families that exist in BOTH Option 1 and Option 2.
    let mut viable_lineages: Vec<(Milestone12LineageKey, Vec<usize>)> = option2_by_lineage
        .into_iter()
        .filter(|(key, _)| option1_lineage_keys.contains_key(key))
        .collect();

    if viable_lineages.is_empty() {
        tracing::warn!(
            "No shared lineage between Option 1 and Option 2 — using best venturi family"
        );
        sort_report_designs(&mut option2_pool);
        option2_pool.truncate(5);
        for (idx, d) in option2_pool.iter_mut().enumerate() {
            d.rank = idx + 1;
        }
    } else {
        // Sort lineage families by their best Option 2 score (descending).
        viable_lineages.sort_by(|a, b| {
            let a_best =
                a.1.iter()
                    .map(|&i| option2_pool[i].score)
                    .fold(f64::NEG_INFINITY, f64::max);
            let b_best =
                b.1.iter()
                    .map(|&i| option2_pool[i].score)
                    .fold(f64::NEG_INFINITY, f64::max);
            b_best.total_cmp(&a_best)
        });

        // Prefer a lineage with >= 5 candidates; fall back to the highest-scoring.
        let winning_indices = viable_lineages
            .iter()
            .find(|(_, indices)| indices.len() >= 5)
            .or(viable_lineages.first())
            .ok_or("No shared Option 1→Option 2 family retained derived venturi candidates")?
            .1
            .clone();

        // Extract the winning lineage's designs from option2_pool by swapping
        // them into a new vec, avoiding clone.
        let mut kept: Vec<bool> = vec![false; option2_pool.len()];
        for &i in &winning_indices {
            kept[i] = true;
        }
        // Filter option2_pool to winning lineage without cloning.
        let mut lineage_pool: Vec<Milestone12ReportDesign> = option2_pool
            .into_iter()
            .enumerate()
            .filter_map(|(i, d)| if kept[i] { Some(d) } else { None })
            .collect();
        sort_report_designs(&mut lineage_pool);
        lineage_pool.truncate(5);
        for (idx, d) in lineage_pool.iter_mut().enumerate() {
            d.rank = idx + 1;
        }
        option2_pool = lineage_pool;
    }

    // Free Option 1 pools — no longer needed after lineage selection.
    drop(option1_pool);
    drop(option1_from_disk);

    let option2_pool_len = option2_pool.len();
    if option2_pool.is_empty() {
        return Err("No eligible Option 2 candidates".into());
    }

    tracing::info!(
        "Option 2 rank-1: {} (score={:.4}  \u{03c3}={:.3})",
        option2_pool[0].topology_short_code(),
        option2_pool[0].score,
        option2_pool[0].metrics.cavitation_number
    );

    // ── Validate ─────────────────────────────────────────────────────────────
    for design in &option2_pool {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option2Derived)?;
    }

    // ── Save track outputs ───────────────────────────────────────────────────
    save_top5_report_json(
        &option2_pool,
        &out_dir.join("two_concept_option2_venturi_top5.json"),
    )?;
    tracing::info!("Saved: two_concept_option2_venturi_top5.json");

    // Save scored seeds for GA warm-start, then drop to free memory.
    {
        scored_seeds.sort_by(|a, b| b.0.total_cmp(&a.0));
        let seed_candidates: Vec<&BlueprintCandidate> =
            scored_seeds.iter().take(150).map(|(_, c)| c).collect();
        save_json_pretty(&seed_candidates, &out_dir.join("option2_scored_seeds.json"))?;
        tracing::info!(
            "Saved: option2_scored_seeds.json ({} seeds)",
            seed_candidates.len()
        );
    }
    drop(scored_seeds);

    // ── Robustness screening ─────────────────────────────────────────────────
    let option2_robustness: Vec<RobustnessReport> = if is_fast {
        tracing::info!("Fast mode: skipping robustness sweep");
        Vec::new()
    } else {
        tracing::info!("Running robustness screening on Option 2 shortlist ...");
        let robustness: Vec<RobustnessReport> = option2_pool
            .iter()
            .map(|d| {
                robustness_sweep_blueprint(
                    &d.candidate,
                    OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
                    &STANDARD_PERTURBATIONS,
                )
            })
            .collect();
        let robust_count = robustness.iter().filter(|r| r.is_robust).count();
        tracing::info!("Robust: {}/{}", robust_count, robustness.len());
        let json = serde_json::to_string_pretty(&robustness)?;
        std::fs::write(out_dir.join("option2_combined_robustness_top5.json"), json)?;
        robustness
    };

    // ── Figure 5 ─────────────────────────────────────────────────────────────
    save_figure(
        option2_pool[0].candidate.blueprint(),
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Figure 5 (Option 2 combined selective venturi)",
    );

    // ── Assemble report ──────────────────────────────────────────────────────
    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let ga_top = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    if !ga_top.is_empty() {
        let canonical = workspace_root.join("report").join("milestone12_results.md");
        write_milestone12_results(
            total_candidates,
            option1_ranked.len(),
            option2_pool_len,
            &option1_ranked,
            &option2_pool,
            &ga_top[0],
            &[],
            &option2_robustness,
            &canonical,
        )?;
        tracing::info!("Updated canonical results: {}", canonical.display());

        let narrative_artifacts = write_milestone12_narrative_report(
            &workspace_root,
            &canonical,
            &Milestone12NarrativeInput {
                total_candidates,
                option1_pool_len: option1_ranked.len(),
                option2_pool_len,
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_pool,
                ga_top: &ga_top,
                option2_pool_all: &option2_pool,
                ga_pool_all: &ga_top,
                validation_rows: &[],
                option2_robustness: &option2_robustness,
                ga_best_per_gen: &[],
                fast_mode: is_fast,
            },
        )?;
        tracing::info!(
            "Narrative report: {}",
            narrative_artifacts.narrative_path.display()
        );
    } else {
        tracing::info!("Skipping full report assembly — run milestone12_ga first for GA track");
    }

    tracing::info!(
        "\n=== Option 2 complete in {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}
