//! Milestone 12 — GA: Blueprint Genetic Refinement.
//!
//! Warm-starts a genetic algorithm from Option 2's top scored seeds and
//! refines the venturi-cavitation design space.  Loads seeds from
//! `option2_scored_seeds.json` (written by `milestone12_option2`); falls back
//! to rebuilding from the candidate space if not available.
//!
//! Produces:
//! - `report/milestone12/ga_hydrosdt_top5.json`
//! - `report/figures/top_hydrosdt_schematic.svg`
//! - Regenerates the narrative report.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example milestone12_ga --no-default-features
//! M12_FAST=1 cargo run -p cfd-optim --example milestone12_ga --no-default-features
//! ```

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, compute_blueprint_report_metrics,
    evaluate_blueprint_candidate, is_milestone12_lineage_topology, load_top5_report_json,
    milestone12_lineage_key, save_blueprint_schematic_svg, save_top5_report_json, score_candidate,
    sort_report_designs, validate_milestone12_candidate, write_milestone12_narrative_report,
    write_milestone12_results, BlueprintCandidate, BlueprintGeneticOptimizer,
    Milestone12LineageKey, Milestone12NarrativeInput, Milestone12ReportDesign, Milestone12Stage,
    OptimMode, OptimizationGoal, SdtWeights,
};
use rayon::prelude::*;
use std::path::{Path, PathBuf};
use std::time::Instant;

// ── Helpers ──────────────────────────────────────────────────────────────────

fn is_selective_report_topology(candidate: &BlueprintCandidate) -> bool {
    is_milestone12_lineage_topology(candidate)
}

fn ga_matches_lineage_sequence(
    candidate: &BlueprintCandidate,
    selected_key: &Milestone12LineageKey,
) -> bool {
    candidate
        .topology_spec()
        .is_ok_and(|spec| spec.stage_sequence_label() == selected_key.stage_sequence_label())
}

fn save_figure(blueprint: &cfd_schematics::NetworkBlueprint, path: &Path, label: &str) {
    match save_blueprint_schematic_svg(blueprint, path) {
        Ok(()) => tracing::info!("      \u{2713} {label}  \u{2192}  {}", path.display()),
        Err(e) => tracing::warn!("      \u{2717} {label}  FAILED: {e}"),
    }
}

fn fast_mode() -> bool {
    std::env::var("M12_FAST")
        .ok()
        .map(|v| {
            let s = v.trim().to_ascii_lowercase();
            s == "1" || s == "true" || s == "yes"
        })
        .unwrap_or(false)
}

fn fast_env(name: &str, default: usize) -> usize {
    std::env::var(name)
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(default)
        .max(1)
}

/// Try to load scored seeds saved by milestone12_option2.
fn load_seeds(out_dir: &Path) -> Result<Vec<BlueprintCandidate>, Box<dyn std::error::Error>> {
    let path = out_dir.join("option2_scored_seeds.json");
    if !path.exists() {
        return Ok(Vec::new());
    }
    let contents = std::fs::read_to_string(&path)?;
    let seeds: Vec<BlueprintCandidate> = serde_json::from_str(&contents)?;
    Ok(seeds)
}

/// Fall back: rebuild seeds from the candidate space by evaluating venturi
/// candidates and taking the top-scoring ones.
fn rebuild_seeds(take: usize) -> Result<Vec<BlueprintCandidate>, Box<dyn std::error::Error>> {
    tracing::info!("Rebuilding GA seeds from candidate space ...");
    let all = build_milestone12_blueprint_candidate_space()?;
    let ga_mode = OptimMode::HydrodynamicCavitationSDT;
    let weights = SdtWeights::default();

    let mut scored: Vec<(f64, BlueprintCandidate)> = all
        .into_par_iter()
        .filter(|c| {
            c.topology_spec()
                .is_ok_and(cfd_schematics::BlueprintTopologySpec::has_venturi)
                && is_selective_report_topology(c)
        })
        .filter_map(|candidate| {
            let _eval = evaluate_blueprint_candidate(&candidate).ok()?;
            let metrics = compute_blueprint_report_metrics(&candidate).ok()?;
            if !metrics.pressure_feasible || !metrics.plate_fits || !metrics.fda_main_compliant {
                return None;
            }
            let score = score_candidate(&metrics, ga_mode, &weights);
            if score <= 0.0 {
                return None;
            }
            Some((score, candidate))
        })
        .collect();

    scored.sort_by(|a, b| b.0.total_cmp(&a.0));
    Ok(scored.into_iter().take(take).map(|(_, c)| c).collect())
}

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let _ = tracing_subscriber::fmt()
        .with_env_filter(
            tracing_subscriber::EnvFilter::builder()
                .with_default_directive(tracing::Level::INFO.into())
                .from_env_lossy(),
        )
        .try_init();

    let workspace_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("cfd-optim crate has a parent")
        .parent()
        .expect("crates/ has a workspace root")
        .to_path_buf();

    let out_dir = workspace_root.join("report").join("milestone12");
    let figures_dir = workspace_root.join("report").join("figures");
    std::fs::create_dir_all(&out_dir)?;
    std::fs::create_dir_all(&figures_dir)?;

    tracing::info!("=== Milestone 12 — GA: Blueprint Genetic Refinement ===\n");
    let overall_start = Instant::now();
    let is_fast = fast_mode();
    let ga_mode = OptimMode::HydrodynamicCavitationSDT;
    let weights = SdtWeights::default();

    let seed_take = if is_fast {
        fast_env("M12_FAST_GA_SEED_TAKE", 24)
    } else {
        150
    };
    let ga_population = if is_fast {
        fast_env("M12_FAST_GA_POPULATION", 12)
    } else {
        80
    };
    let ga_generations = if is_fast {
        fast_env("M12_FAST_GA_GENERATIONS", 8)
    } else {
        120
    };

    // ── Load or rebuild seeds ────────────────────────────────────────────────
    let mut seeds = load_seeds(&out_dir)?;
    if seeds.is_empty() {
        seeds = rebuild_seeds(seed_take)?;
    }
    seeds.truncate(seed_take);

    if seeds.is_empty() {
        return Err(
            "GA received no seeds — run milestone12_option2 first or ensure candidates are feasible"
                .into(),
        );
    }
    tracing::info!("GA warm-start: {} seeds", seeds.len());

    // ── Determine lineage key from Option 2 results ──────────────────────────
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;
    let selected_lineage_key: Option<Milestone12LineageKey> = option2_ranked
        .first()
        .and_then(|d| milestone12_lineage_key(&d.candidate));

    // ── Run GA ───────────────────────────────────────────────────────────────
    tracing::info!(
        "Running HydroSDT GA (pop={}, gen={}) ...",
        ga_population,
        ga_generations
    );
    let ga_start = Instant::now();

    let ga_result = BlueprintGeneticOptimizer::new(OptimizationGoal::BlueprintGeneticRefinement)
        .with_seeds(seeds.clone())
        .with_population(ga_population)
        .with_max_generations(ga_generations)
        .with_top_k(5)
        .run()?;

    tracing::info!(
        "GA complete in {:.1}s ({} total candidates evaluated)",
        ga_start.elapsed().as_secs_f64(),
        ga_result.all_candidates.len()
    );

    // ── Build GA report pool ─────────────────────────────────────────────────
    let mut ga_report_top: Vec<Milestone12ReportDesign> = ga_result
        .top_candidates
        .iter()
        .filter(|ranked| {
            if let Some(ref key) = selected_lineage_key {
                is_selective_report_topology(&ranked.candidate)
                    && ga_matches_lineage_sequence(&ranked.candidate, key)
            } else {
                is_selective_report_topology(&ranked.candidate)
            }
        })
        .filter_map(|ranked| {
            Milestone12ReportDesign::from_blueprint_candidate(
                ranked.rank,
                ranked.candidate.clone(),
                ranked.evaluation.score,
            )
            .ok()
        })
        .collect();

    // Fallback: if lineage filter is too strict, use seeds directly.
    if ga_report_top.is_empty() {
        tracing::warn!("GA top_candidates produced no lineage-matched designs; falling back to seeds");
        ga_report_top = seeds
            .iter()
            .filter(|c| {
                is_selective_report_topology(c)
                    && selected_lineage_key.as_ref().map_or(true, |key| {
                        ga_matches_lineage_sequence(c, key)
                    })
            })
            .filter_map(|candidate| {
                let metrics = compute_blueprint_report_metrics(candidate).ok()?;
                let score = score_candidate(&metrics, ga_mode, &weights);
                Some(Milestone12ReportDesign::new(
                    0,
                    candidate.clone(),
                    metrics,
                    score,
                ))
            })
            .collect();
        sort_report_designs(&mut ga_report_top);
        ga_report_top.truncate(5);
    }

    // Second fallback: use Option 2 designs directly.
    if ga_report_top.is_empty() && !option2_ranked.is_empty() {
        tracing::warn!("Using Option 2 designs as GA fallback");
        ga_report_top = option2_ranked
            .iter()
            .take(5)
            .map(|d| {
                Milestone12ReportDesign::new(
                    0,
                    d.candidate.clone(),
                    d.metrics.clone(),
                    score_candidate(&d.metrics, ga_mode, &weights),
                )
            })
            .collect();
        sort_report_designs(&mut ga_report_top);
    }

    if ga_report_top.is_empty() {
        return Err("GA produced no viable designs".into());
    }

    for (idx, d) in ga_report_top.iter_mut().enumerate() {
        d.rank = idx + 1;
    }

    tracing::info!(
        "GA rank-1: {} (score={:.4}  \u{03c3}={:.3}  cancer_cav={:.3})",
        ga_report_top[0].candidate.id,
        ga_report_top[0].score,
        ga_report_top[0].metrics.cavitation_number,
        ga_report_top[0].metrics.cancer_targeted_cavitation
    );

    // ── Validate ─────────────────────────────────────────────────────────────
    for design in &ga_report_top {
        validate_milestone12_candidate(&design.candidate, Milestone12Stage::GaRefined)?;
    }

    // ── Save track outputs ───────────────────────────────────────────────────
    save_top5_report_json(&ga_report_top, &out_dir.join("ga_hydrosdt_top5.json"))?;
    tracing::info!("Saved: ga_hydrosdt_top5.json");

    // ── Figure 6 ─────────────────────────────────────────────────────────────
    save_figure(
        ga_report_top[0].candidate.blueprint(),
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Figure 6 (HydroSDT GA rank-1)",
    );

    // ── Assemble report ──────────────────────────────────────────────────────
    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;

    // Build the full GA pool for the narrative (all candidates from GA run).
    let ga_pool_all: Vec<Milestone12ReportDesign> = ga_result
        .all_candidates
        .iter()
        .filter(|ranked| is_selective_report_topology(&ranked.candidate))
        .filter_map(|ranked| {
            Milestone12ReportDesign::from_blueprint_candidate(
                ranked.rank,
                ranked.candidate.clone(),
                ranked.evaluation.score,
            )
            .ok()
        })
        .collect();

    if !option2_ranked.is_empty() {
        let canonical = workspace_root.join("report").join("milestone12_results.md");
        write_milestone12_results(
            0,
            option1_ranked.len(),
            option2_ranked.len(),
            &option1_ranked,
            &option2_ranked,
            &ga_report_top[0],
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
                option1_pool_len: option1_ranked.len(),
                option2_pool_len: option2_ranked.len(),
                option1_ranked: &option1_ranked,
                option2_ranked: &option2_ranked,
                ga_top: &ga_report_top,
                option2_pool_all: &option2_ranked,
                ga_pool_all: &ga_pool_all,
                validation_rows: &[],
                option2_robustness: &[],
                ga_best_per_gen: &ga_result.best_per_generation,
                fast_mode: is_fast,
            },
        )?;
        tracing::info!(
            "Narrative report: {}",
            narrative_artifacts.narrative_path.display()
        );
    } else {
        tracing::info!("Skipping report assembly — run milestone12_option2 first");
    }

    tracing::info!(
        "\n=== GA complete in {:.1}s ===",
        overall_start.elapsed().as_secs_f64()
    );

    Ok(())
}
