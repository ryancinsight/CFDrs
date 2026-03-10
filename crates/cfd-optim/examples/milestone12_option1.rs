//! Milestone 12 — Option 1: Selective Acoustic Residence Separation.
//!
//! Evaluates non-venturi selective-routing topologies using the blueprint-native
//! `EvaluatedPool` scoring system:
//!
//!   score = residence × separation × safety  (multiplicative)
//!
//! This correctly ranks Tri→Tri→Tri above Tri→Tri when its additional
//! Zweifach–Fung splitting stage achieves higher cancer-center enrichment
//! (~0.75+ vs ~0.556), because the multiplicative product amplifies the
//! separation advantage rather than diluting it across additive weight budgets.
//!
//! Produces:
//! - `report/milestone12/two_concept_option1_ultrasound_top5.json`
//! - `report/figures/selected_option1_schematic.svg`
//! - Regenerates the narrative report (using Option 2 / GA data from prior runs
//!   if available on disk).
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_option1 --no-default-features
//! # fast iteration
//! M12_FAST=1 cargo run -p cfd-optim --example milestone12_option1 --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, fast_env, fast_mode, init_tracing,
    is_selective_report_topology, load_top5_report_json, resolve_output_directories, save_figure,
    save_json_pretty, save_top5_report_json, validate_milestone12_candidate,
    write_milestone12_narrative_report, write_milestone12_results, BlueprintCandidate, EvaluatedPool,
    Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage, OptimizationGoal,
};
use serde::Serialize;
use std::time::Instant;

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

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    init_tracing();
    let (workspace_root, out_dir, figures_dir) = resolve_output_directories()?;

    tracing::info!(
        "=== Milestone 12 — Option 1: Selective Acoustic (blueprint-native scoring) ===\n"
    );
    let overall_start = Instant::now();
    let is_fast = fast_mode();

    // ── Build candidate space ────────────────────────────────────────────────
    let t_build = Instant::now();
    let raw_candidates = build_milestone12_blueprint_candidate_space()?;
    tracing::info!(
        "Candidate space: {} blueprints ({:.1}s)",
        raw_candidates.len(),
        t_build.elapsed().as_secs_f64()
    );

    // ── Filter to non-venturi selective topologies ───────────────────────────
    let acoustic_candidates: Vec<BlueprintCandidate> = raw_candidates
        .into_iter()
        .filter(|c| {
            !c.topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
                && is_selective_report_topology(c)
        })
        .collect();

    let eval_cap = if is_fast {
        fast_env("M12_FAST_ACOUSTIC_EVAL_MAX", 200)
    } else {
        acoustic_candidates.len()
    };
    let to_evaluate: Vec<BlueprintCandidate> =
        acoustic_candidates.into_iter().take(eval_cap).collect();
    tracing::info!(
        "Non-venturi selective candidates: {} (evaluating via EvaluatedPool)",
        to_evaluate.len()
    );

    // ── Evaluate and rank using blueprint-native multiplicative scoring ──────
    //
    // score = residence × separation × safety
    //   residence  = treatment_residence_time_s × treatment_flow_fraction
    //   separation = separation_efficiency × cancer_center_fraction
    //   safety     = max(0, main_channel_margin)
    //
    // This correctly ranks Tri→Tri→Tri above Tri→Tri because the additional
    // Zweifach–Fung stage yields ~0.75+ cancer_center_fraction vs ~0.556.
    let t_eval = Instant::now();
    let pool = EvaluatedPool::from_candidates(&to_evaluate);
    tracing::info!(
        "Pool: {} / {} evaluated in {:.1}s ({:.1} MB heap)",
        pool.len(),
        to_evaluate.len(),
        t_eval.elapsed().as_secs_f64(),
        pool.heap_bytes() as f64 / (1024.0 * 1024.0),
    );

    // Count eligible entries without materializing them (avoids OOM).
    let shortlist_size = if is_fast { 3 } else { 5 };
    let goal = OptimizationGoal::AsymmetricSplitResidenceSeparation;
    let option1_pool_len = pool.count_eligible(goal);
    tracing::info!("Option 1 eligible (score > 0): {}", option1_pool_len);

    // Only materialise the shortlist — NOT the full eligible pool.
    let top_results = pool.top_k(shortlist_size, goal);

    // ── Convert top results to report designs ────────────────────────────────
    //
    // Milestone12ReportDesign still requires SdtMetrics for the downstream
    // narrative report.  We compute legacy metrics only for the shortlisted
    // results — the expensive ranking already happened via the pool.
    let mut option1_ranked: Vec<Milestone12ReportDesign> = match top_results {
        Ok(results) => results
            .into_iter()
            .enumerate()
            .filter_map(|(idx, eval)| {
                let candidate = to_evaluate
                    .iter()
                    .find(|c| c.id == eval.candidate_id)?
                    .clone();
                match Milestone12ReportDesign::from_blueprint_candidate(
                    idx + 1,
                    candidate,
                    eval.score,
                ) {
                    Ok(design) => Some(design),
                    Err(e) => {
                        tracing::warn!(
                            "Failed to compute report metrics for {}: {e}",
                            eval.candidate_id
                        );
                        None
                    }
                }
            })
            .collect(),
        Err(e) => {
            tracing::warn!("No eligible Option 1 candidates: {e}");
            Vec::new()
        }
    };

    for (idx, design) in option1_ranked.iter_mut().enumerate() {
        design.rank = idx + 1;
    }

    if let Some(best) = option1_ranked.first() {
        tracing::info!(
            "Option 1 rank-1: {} (score={:.4})",
            best.topology_short_code(),
            best.score
        );
    } else {
        tracing::warn!("No eligible Option 1 candidates under current physics regime");
    }

    // ── Validate ─────────────────────────────────────────────────────────────
    for design in &option1_ranked {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::Option1Base)?;
    }

    // ── Save track outputs ───────────────────────────────────────────────────
    save_top5_report_json(
        &option1_ranked,
        &out_dir.join("two_concept_option1_ultrasound_top5.json"),
    )?;
    tracing::info!("Saved: two_concept_option1_ultrasound_top5.json");

    if let Some(topline) = option1_ranked.first().map(|d| {
        concept_topline(
            "Option 1: Selective acoustic center treatment",
            "AsymmetricSplitResidenceSeparation",
            d,
        )
    }) {
        save_json_pretty(
            &vec![topline],
            &out_dir.join("option1_concept_topline.json"),
        )?;
    }

    // ── Figure 4 ─────────────────────────────────────────────────────────────
    if let Some(best) = option1_ranked.first() {
        save_figure(
            best.candidate.blueprint(),
            &figures_dir.join("selected_option1_schematic.svg"),
            "Figure 4 (Option 1 selective acoustic)",
        );
    }

    // ── Assemble report (load other tracks if available) ─────────────────────
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let ga_top = load_top5_report_json(&out_dir.join("ga_hydrosdt_top5.json"))?;

    if !option2_ranked.is_empty() && !ga_top.is_empty() {
        let canonical = workspace_root.join("report").join("milestone12_results.md");
        write_milestone12_results(
            0, // total_candidates not tracked in this partial run
            option1_pool_len,
            option2_ranked.len(),
            &option1_ranked,
            &option2_ranked,
            &ga_top[0],
            &[],
            &[],
            &canonical,
        )?;
        tracing::info!("Updated canonical results: {}", canonical.display());

        let narrative_artifacts = write_milestone12_narrative_report(
            &workspace_root,
            &canonical,
            &Milestone12NarrativeInput {
                total_candidates: 0,
                option1_pool_len,
                option2_pool_len: option2_ranked.len(),
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_ranked,
                ga_top: &ga_top,
                option2_pool_all: &option2_ranked,
                ga_pool_all: &ga_top,
                validation_rows: &[],
                option2_robustness: &[],
                ga_best_per_gen: &[],
                fast_mode: is_fast,
            },
        )?;
        tracing::info!(
            "Narrative report: {}",
            narrative_artifacts.narrative_path.display()
        );
    } else {
        tracing::info!(
            "Skipping report assembly — run milestone12_option2 and milestone12_ga first"
        );
    }

    tracing::info!(
        "\n=== Option 1 complete in {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}
