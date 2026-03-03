//! Exhaustive scenario-space audit for `cfd-schematics` + `cfd-optim`.
//!
//! This example does two things in one run:
//! 1. Catalogues every preset scenario currently exported by
//!    `cfd_schematics::interface::presets`.
//! 2. Runs a full parametric simulation sweep (`build_candidate_space`) and
//!    reports per-topology/per-mode feasibility and ranking diagnostics.
//!
//! # Run
//! ```bash
//! cargo run -p cfd-optim --example sdt_scenario_audit --no-default-features
//! ```

use cfd_optim::{
    build_candidate_space, compute_metrics, score_candidate, DesignCandidate, OptimMode,
    SdtMetrics, SdtWeights,
};
use cfd_schematics::interface::presets::{
    asymmetric_bifurcation_serpentine_rect, bifurcation_rect, bifurcation_serpentine_rect,
    bifurcation_trifurcation_venturi_rect, bifurcation_venturi_rect,
    cascade_center_trifurcation_rect, cell_separation_rect, constriction_expansion_array_rect,
    double_bifurcation_venturi_rect, double_trifurcation_venturi_rect,
    incremental_filtration_tri_bi_rect, parallel_microchannel_array_rect,
    quad_trifurcation_venturi_rect, serial_double_venturi_rect, serpentine_chain, serpentine_rect,
    spiral_channel_rect, symmetric_bifurcation, symmetric_trifurcation,
    trifurcation_bifurcation_bifurcation_venturi_rect, trifurcation_bifurcation_venturi_rect,
    trifurcation_rect, trifurcation_serpentine_rect, trifurcation_venturi_rect,
    triple_bifurcation_venturi_rect, triple_trifurcation_venturi_rect, venturi_chain, venturi_rect,
    venturi_serpentine_rect,
};
use cfd_schematics::NetworkBlueprint;
use serde::Serialize;
use std::collections::{BTreeMap, HashMap};
use std::path::Path;

const MODE_COUNT: usize = 8;
const TOP_K_PER_MODE: usize = 20;

#[derive(Clone, Copy)]
struct ModeDef {
    key: &'static str,
    label: &'static str,
    mode: OptimMode,
}

#[derive(Clone, Copy)]
struct PresetDef {
    preset_id: &'static str,
    preset_fn: &'static str,
    family: &'static str,
    support: &'static str,
    mapped_to: &'static str,
    build: fn(&str) -> NetworkBlueprint,
}

#[derive(Debug, Clone, Serialize)]
struct PresetRow {
    preset_id: String,
    preset_fn: String,
    family: String,
    support: String,
    mapped_to: String,
    nodes: usize,
    inlets: usize,
    outlets: usize,
    junctions: usize,
    channels: usize,
    total_length_mm: f64,
    venturi_channels: usize,
    serpentine_like_channels: usize,
}

#[derive(Debug, Clone)]
struct EvaluatedRecord {
    candidate: DesignCandidate,
    metrics: SdtMetrics,
    scores: [f64; MODE_COUNT],
}

#[derive(Debug, Clone, Serialize)]
struct ModeSummary {
    mode_key: String,
    mode_label: String,
    feasible_candidates: usize,
    feasible_rate_pct: f64,
    best_score: f64,
    mean_feasible_score: f64,
}

#[derive(Debug, Clone, Serialize)]
struct TopologyModeRow {
    mode_key: String,
    mode_label: String,
    topology_short: String,
    topology_name: String,
    total_candidates: usize,
    evaluated_candidates: usize,
    failed_metrics: usize,
    feasible_count: usize,
    feasible_rate_pct: f64,
    hydro_cav_ready_count: usize,
    best_score: f64,
    median_feasible_score: f64,
    p90_feasible_score: f64,
    mean_feasible_score: f64,
    best_candidate_id: String,
    best_sigma: f64,
    best_hi_per_pass: f64,
    best_bulk_hi_per_pass: f64,
    best_sep3_eff: f64,
    best_wbc_recovery: f64,
    best_local_hct_venturi: f64,
    best_cancer_dose: f64,
    best_venturi_flow_frac: f64,
    best_rbc_venturi_exposure: f64,
    best_delta_p_kpa: f64,
}

#[derive(Debug, Clone, Serialize)]
struct TopCandidateRow {
    mode_key: String,
    mode_label: String,
    rank: usize,
    candidate_id: String,
    topology_short: String,
    topology_name: String,
    score: f64,
    sigma: f64,
    hi_per_pass: f64,
    bulk_hi_per_pass: f64,
    sep3_eff: f64,
    wbc_recovery: f64,
    coverage_pct: f64,
    delta_p_kpa: f64,
    local_hct_venturi: f64,
    cancer_dose_fraction: f64,
    venturi_flow_fraction: f64,
    rbc_venturi_exposure_fraction: f64,
    fda_main: bool,
    pressure_feasible: bool,
    hydro_cav_ready: bool,
}

#[derive(Debug, Serialize)]
struct ScenarioAudit {
    total_presets: usize,
    total_candidates: usize,
    evaluated_candidates: usize,
    failed_metrics: usize,
    mode_summaries: Vec<ModeSummary>,
    preset_catalog: Vec<PresetRow>,
    topology_mode_stats: Vec<TopologyModeRow>,
    top_candidates_by_mode: BTreeMap<String, Vec<TopCandidateRow>>,
}

fn main() {
    let out_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("outputs")
        .join("sdt_scenario_audit");
    std::fs::create_dir_all(&out_dir).expect("failed to create output directory");

    println!("{}", "=".repeat(112));
    println!(
        "  cfd-optim  |  Scenario Audit  |  cfd-schematics surface + exhaustive simulation matrix"
    );
    println!("{}", "=".repeat(112));

    // ── Part 1: Preset-surface audit ─────────────────────────────────────────
    let preset_defs = preset_defs();
    let preset_rows: Vec<PresetRow> = preset_defs.iter().map(inspect_preset).collect();

    let mapped_count = preset_rows
        .iter()
        .filter(|r| !r.mapped_to.is_empty())
        .count();
    let schematics_only_count = preset_rows.len().saturating_sub(mapped_count);

    println!("\n[1/4] Preset-surface audit");
    println!(
        "  - Presets exported by cfd-schematics: {}",
        preset_rows.len()
    );
    println!(
        "  - Mapped into cfd-optim design topologies: {}",
        mapped_count
    );
    println!(
        "  - Schematics-only (not in current optimizer sweep): {}",
        schematics_only_count
    );

    // ── Part 2: Full simulation sweep ────────────────────────────────────────
    println!("\n[2/4] Exhaustive simulation sweep over build_candidate_space()");
    let modes = modes();
    let weights = SdtWeights::default();
    let candidates = build_candidate_space();
    let total_candidates = candidates.len();

    let mut topology_totals: HashMap<String, usize> = HashMap::new();
    let mut topology_failures: HashMap<String, usize> = HashMap::new();
    let mut topology_names: HashMap<String, String> = HashMap::new();
    let mut evaluated = Vec::with_capacity(total_candidates);

    for candidate in candidates {
        let topo_short = candidate.topology.short().to_string();
        topology_names
            .entry(topo_short.clone())
            .or_insert_with(|| candidate.topology.name().to_string());
        *topology_totals.entry(topo_short.clone()).or_insert(0) += 1;

        match compute_metrics(&candidate) {
            Ok(metrics) => {
                let mut scores = [0.0; MODE_COUNT];
                for (i, mode_def) in modes.iter().enumerate() {
                    scores[i] = score_candidate(&metrics, mode_def.mode, &weights);
                }
                evaluated.push(EvaluatedRecord {
                    candidate,
                    metrics,
                    scores,
                });
            }
            Err(_) => {
                *topology_failures.entry(topo_short).or_insert(0) += 1;
            }
        }
    }

    let failed_metrics = total_candidates.saturating_sub(evaluated.len());
    println!("  - Total candidates generated: {total_candidates}");
    println!("  - Successfully evaluated: {}", evaluated.len());
    println!("  - Metrics evaluation failures: {failed_metrics}");

    // Index evaluated records by topology short code for fast aggregation.
    let mut by_topology: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, rec) in evaluated.iter().enumerate() {
        by_topology
            .entry(rec.candidate.topology.short().to_string())
            .or_default()
            .push(idx);
    }

    // ── Part 3: Per-mode aggregation ─────────────────────────────────────────
    println!("\n[3/4] Aggregating mode-wise feasibility and ranking diagnostics");
    let mut mode_summaries = Vec::with_capacity(MODE_COUNT);
    let mut topology_mode_stats = Vec::new();
    let mut top_candidates_by_mode: BTreeMap<String, Vec<TopCandidateRow>> = BTreeMap::new();

    let mut topology_keys: Vec<String> = topology_totals.keys().cloned().collect();
    topology_keys.sort();

    for (mode_idx, mode_def) in modes.iter().enumerate() {
        let mut feasible_total = 0usize;
        let mut feasible_scores_global = Vec::new();

        for topo_short in &topology_keys {
            let total = *topology_totals.get(topo_short).unwrap_or(&0);
            let failed = *topology_failures.get(topo_short).unwrap_or(&0);
            let rec_ids: &[usize] = by_topology
                .get(topo_short)
                .map(std::vec::Vec::as_slice)
                .unwrap_or(&[]);
            let evaluated_count = rec_ids.len();

            let mut feasible_scores = Vec::new();
            let mut hydro_cav_ready_count = 0usize;
            let mut best_idx: Option<usize> = None;
            let mut best_score = 0.0_f64;

            for &rec_idx in rec_ids {
                let rec = &evaluated[rec_idx];
                let score = rec.scores[mode_idx];
                if score > 0.0 {
                    feasible_scores.push(score);
                    feasible_scores_global.push(score);
                    feasible_total += 1;
                    if is_hydrodynamic_cavitation_ready(&rec.metrics) {
                        hydro_cav_ready_count += 1;
                    }
                    if score > best_score {
                        best_score = score;
                        best_idx = Some(rec_idx);
                    }
                }
            }

            feasible_scores.sort_by(|a, b| a.total_cmp(b));
            let feasible_count = feasible_scores.len();
            let feasible_rate = if total > 0 {
                feasible_count as f64 / total as f64 * 100.0
            } else {
                0.0
            };
            let median = quantile(&feasible_scores, 0.50);
            let p90 = quantile(&feasible_scores, 0.90);
            let mean = if feasible_scores.is_empty() {
                0.0
            } else {
                feasible_scores.iter().sum::<f64>() / feasible_scores.len() as f64
            };

            let topology_name = topology_names
                .get(topo_short)
                .cloned()
                .unwrap_or_else(|| "Unknown".to_string());

            if let Some(idx) = best_idx {
                let best_rec = &evaluated[idx];
                topology_mode_stats.push(TopologyModeRow {
                    mode_key: mode_def.key.to_string(),
                    mode_label: mode_def.label.to_string(),
                    topology_short: topo_short.clone(),
                    topology_name,
                    total_candidates: total,
                    evaluated_candidates: evaluated_count,
                    failed_metrics: failed,
                    feasible_count,
                    feasible_rate_pct: feasible_rate,
                    hydro_cav_ready_count,
                    best_score,
                    median_feasible_score: median,
                    p90_feasible_score: p90,
                    mean_feasible_score: mean,
                    best_candidate_id: best_rec.candidate.id.clone(),
                    best_sigma: best_rec.metrics.cavitation_number,
                    best_hi_per_pass: best_rec.metrics.hemolysis_index_per_pass,
                    best_bulk_hi_per_pass: best_rec.metrics.bulk_hemolysis_index_per_pass,
                    best_sep3_eff: best_rec.metrics.three_pop_sep_efficiency,
                    best_wbc_recovery: best_rec.metrics.wbc_recovery,
                    best_local_hct_venturi: best_rec.metrics.local_hematocrit_venturi,
                    best_cancer_dose: best_rec.metrics.cancer_dose_fraction,
                    best_venturi_flow_frac: best_rec.metrics.venturi_flow_fraction,
                    best_rbc_venturi_exposure: best_rec.metrics.rbc_venturi_exposure_fraction,
                    best_delta_p_kpa: best_rec.metrics.total_pressure_drop_pa * 1e-3,
                });
            } else {
                topology_mode_stats.push(TopologyModeRow {
                    mode_key: mode_def.key.to_string(),
                    mode_label: mode_def.label.to_string(),
                    topology_short: topo_short.clone(),
                    topology_name,
                    total_candidates: total,
                    evaluated_candidates: evaluated_count,
                    failed_metrics: failed,
                    feasible_count,
                    feasible_rate_pct: feasible_rate,
                    hydro_cav_ready_count,
                    best_score: 0.0,
                    median_feasible_score: 0.0,
                    p90_feasible_score: 0.0,
                    mean_feasible_score: 0.0,
                    best_candidate_id: String::new(),
                    best_sigma: f64::INFINITY,
                    best_hi_per_pass: 0.0,
                    best_bulk_hi_per_pass: 0.0,
                    best_sep3_eff: 0.0,
                    best_wbc_recovery: 0.0,
                    best_local_hct_venturi: 0.0,
                    best_cancer_dose: 0.0,
                    best_venturi_flow_frac: 0.0,
                    best_rbc_venturi_exposure: 0.0,
                    best_delta_p_kpa: 0.0,
                });
            }
        }

        let mut ranked: Vec<(usize, f64)> = evaluated
            .iter()
            .enumerate()
            .filter_map(|(idx, rec)| {
                let score = rec.scores[mode_idx];
                (score > 0.0).then_some((idx, score))
            })
            .collect();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1));

        let top_rows: Vec<TopCandidateRow> = ranked
            .iter()
            .take(TOP_K_PER_MODE)
            .enumerate()
            .map(|(i, (idx, score))| {
                let rec = &evaluated[*idx];
                TopCandidateRow {
                    mode_key: mode_def.key.to_string(),
                    mode_label: mode_def.label.to_string(),
                    rank: i + 1,
                    candidate_id: rec.candidate.id.clone(),
                    topology_short: rec.candidate.topology.short().to_string(),
                    topology_name: rec.candidate.topology.name().to_string(),
                    score: *score,
                    sigma: rec.metrics.cavitation_number,
                    hi_per_pass: rec.metrics.hemolysis_index_per_pass,
                    bulk_hi_per_pass: rec.metrics.bulk_hemolysis_index_per_pass,
                    sep3_eff: rec.metrics.three_pop_sep_efficiency,
                    wbc_recovery: rec.metrics.wbc_recovery,
                    coverage_pct: rec.metrics.well_coverage_fraction * 100.0,
                    delta_p_kpa: rec.metrics.total_pressure_drop_pa * 1e-3,
                    local_hct_venturi: rec.metrics.local_hematocrit_venturi,
                    cancer_dose_fraction: rec.metrics.cancer_dose_fraction,
                    venturi_flow_fraction: rec.metrics.venturi_flow_fraction,
                    rbc_venturi_exposure_fraction: rec.metrics.rbc_venturi_exposure_fraction,
                    fda_main: rec.metrics.fda_main_compliant,
                    pressure_feasible: rec.metrics.pressure_feasible,
                    hydro_cav_ready: is_hydrodynamic_cavitation_ready(&rec.metrics),
                }
            })
            .collect();
        top_candidates_by_mode.insert(mode_def.key.to_string(), top_rows);

        let feasible_rate_pct = if total_candidates > 0 {
            feasible_total as f64 / total_candidates as f64 * 100.0
        } else {
            0.0
        };
        let best_score = ranked.first().map(|(_, s)| *s).unwrap_or(0.0);
        let mean_feasible_score = if feasible_scores_global.is_empty() {
            0.0
        } else {
            feasible_scores_global.iter().sum::<f64>() / feasible_scores_global.len() as f64
        };
        mode_summaries.push(ModeSummary {
            mode_key: mode_def.key.to_string(),
            mode_label: mode_def.label.to_string(),
            feasible_candidates: feasible_total,
            feasible_rate_pct,
            best_score,
            mean_feasible_score,
        });
    }

    // ── Part 4: Persist audit outputs ────────────────────────────────────────
    println!("\n[4/4] Writing audit artifacts");
    let audit = ScenarioAudit {
        total_presets: preset_rows.len(),
        total_candidates,
        evaluated_candidates: evaluated.len(),
        failed_metrics,
        mode_summaries: mode_summaries.clone(),
        preset_catalog: preset_rows.clone(),
        topology_mode_stats: topology_mode_stats.clone(),
        top_candidates_by_mode: top_candidates_by_mode.clone(),
    };

    let json_path = out_dir.join("scenario_audit.json");
    let md_path = out_dir.join("scenario_audit.md");
    let preset_csv_path = out_dir.join("preset_catalog.csv");
    let topology_csv_path = out_dir.join("topology_mode_stats.csv");

    std::fs::write(
        &json_path,
        serde_json::to_string_pretty(&audit).expect("failed to serialise audit json"),
    )
    .expect("failed to write audit json");

    write_preset_catalog_csv(&preset_rows, &preset_csv_path)
        .expect("failed to write preset catalog csv");
    write_topology_mode_csv(&topology_mode_stats, &topology_csv_path)
        .expect("failed to write topology stats csv");

    for mode in modes {
        let mode_csv = out_dir.join(format!("top_candidates_{}.csv", mode.key));
        let rows = top_candidates_by_mode
            .get(mode.key)
            .cloned()
            .unwrap_or_default();
        write_top_candidates_csv(&rows, &mode_csv).expect("failed to write top candidates csv");
    }

    write_markdown_report(
        &md_path,
        &preset_rows,
        &mode_summaries,
        &topology_mode_stats,
        &top_candidates_by_mode,
    )
    .expect("failed to write markdown report");

    let combined_mode_idx = modes
        .iter()
        .position(|m| m.key == "combined_sdt_leuka")
        .unwrap_or(MODE_COUNT.saturating_sub(1));
    let mut combined_ranked: Vec<(f64, &EvaluatedRecord)> = evaluated
        .iter()
        .filter_map(|rec| {
            let score = rec.scores[combined_mode_idx];
            (score > 0.0).then_some((score, rec))
        })
        .collect();
    combined_ranked.sort_by(|a, b| b.0.total_cmp(&a.0));

    let cif_cct_md_path = out_dir.join("cif_cct_stepwise.md");
    write_cif_cct_stepwise_markdown(&cif_cct_md_path, &combined_ranked)
        .expect("failed to write cif/cct stepwise markdown");

    let report_path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .and_then(|p| p.parent())
        .map(|p| p.join("report").join("milestone12_scenario_audit.md"));
    if let Some(path) = report_path {
        write_cif_cct_stepwise_markdown(&path, &combined_ranked)
            .expect("failed to write report milestone12 scenario markdown");
    }

    println!("  Saved: {}", json_path.display());
    println!("  Saved: {}", md_path.display());
    println!("  Saved: {}", cif_cct_md_path.display());
    println!("  Saved: {}", preset_csv_path.display());
    println!("  Saved: {}", topology_csv_path.display());
    for mode in modes {
        println!(
            "  Saved: {}",
            out_dir
                .join(format!("top_candidates_{}.csv", mode.key))
                .display()
        );
    }
    println!("\n{}", "=".repeat(112));
    println!("  Completed scenario audit.");
    println!(
        "  Note: `AdaptiveTree` is GA-only and generated via `DesignCandidate::to_blueprint`."
    );
    println!("{}", "=".repeat(112));
}

fn modes() -> [ModeDef; MODE_COUNT] {
    [
        ModeDef {
            key: "sdt_therapy",
            label: "SDT Therapy (Selective Sep + HI + Cav)",
            mode: OptimMode::SdtTherapy,
        },
        ModeDef {
            key: "sdt_cavitation",
            label: "SDT Cavitation",
            mode: OptimMode::SdtCavitation,
        },
        ModeDef {
            key: "three_pop",
            label: "Three-Pop Separation",
            mode: OptimMode::ThreePopSeparation,
        },
        ModeDef {
            key: "cell_separation",
            label: "Cell Separation",
            mode: OptimMode::CellSeparation,
        },
        ModeDef {
            key: "uniform_exposure",
            label: "Uniform Exposure",
            mode: OptimMode::UniformExposure,
        },
        ModeDef {
            key: "pediatric_3kg",
            label: "Pediatric Leukapheresis (3 kg)",
            mode: OptimMode::PediatricLeukapheresis {
                patient_weight_kg: 3.0,
            },
        },
        ModeDef {
            key: "hydro_cav_sdt",
            label: "Hydrodynamic Cavitation SDT",
            mode: OptimMode::HydrodynamicCavitationSDT,
        },
        ModeDef {
            key: "combined_sdt_leuka",
            label: "Combined SDT + Leukapheresis (Milestone 12)",
            mode: OptimMode::CombinedSdtLeukapheresis {
                leuka_weight: 0.5,
                sdt_weight: 0.5,
                patient_weight_kg: 3.0,
            },
        },
    ]
}

fn preset_defs() -> Vec<PresetDef> {
    vec![
        PresetDef {
            preset_id: "venturi_chain",
            preset_fn: "venturi_chain",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_venturi_chain,
        },
        PresetDef {
            preset_id: "venturi_rect",
            preset_fn: "venturi_rect",
            family: "base_rectangular",
            support: "parametric+ga",
            mapped_to: "SV",
            build: preset_venturi_rect,
        },
        PresetDef {
            preset_id: "symmetric_bifurcation",
            preset_fn: "symmetric_bifurcation",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_symmetric_bifurcation,
        },
        PresetDef {
            preset_id: "bifurcation_rect",
            preset_fn: "bifurcation_rect",
            family: "base_rectangular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_bifurcation_rect,
        },
        PresetDef {
            preset_id: "symmetric_trifurcation",
            preset_fn: "symmetric_trifurcation",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_symmetric_trifurcation,
        },
        PresetDef {
            preset_id: "trifurcation_rect",
            preset_fn: "trifurcation_rect",
            family: "base_rectangular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_trifurcation_rect,
        },
        PresetDef {
            preset_id: "serpentine_chain",
            preset_fn: "serpentine_chain",
            family: "base_circular",
            support: "schematics-only",
            mapped_to: "",
            build: preset_serpentine_chain,
        },
        PresetDef {
            preset_id: "serpentine_rect",
            preset_fn: "serpentine_rect",
            family: "base_rectangular",
            support: "parametric+ga",
            mapped_to: "SG",
            build: preset_serpentine_rect,
        },
        PresetDef {
            preset_id: "venturi_serpentine_rect",
            preset_fn: "venturi_serpentine_rect",
            family: "series_composite",
            support: "parametric+ga",
            mapped_to: "VS",
            build: preset_venturi_serpentine_rect,
        },
        PresetDef {
            preset_id: "serial_double_venturi_rect",
            preset_fn: "serial_double_venturi_rect",
            family: "series_composite",
            support: "parametric+ga",
            mapped_to: "S2",
            build: preset_serial_double_venturi_rect,
        },
        PresetDef {
            preset_id: "bifurcation_venturi_rect",
            preset_fn: "bifurcation_venturi_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "BV",
            build: preset_bifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "bifurcation_serpentine_rect",
            preset_fn: "bifurcation_serpentine_rect",
            family: "bifurcation_composite",
            support: "parametric+ga",
            mapped_to: "BS",
            build: preset_bifurcation_serpentine_rect,
        },
        PresetDef {
            preset_id: "trifurcation_venturi_rect",
            preset_fn: "trifurcation_venturi_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "TV",
            build: preset_trifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "trifurcation_serpentine_rect",
            preset_fn: "trifurcation_serpentine_rect",
            family: "trifurcation_composite",
            support: "parametric+ga",
            mapped_to: "TS",
            build: preset_trifurcation_serpentine_rect,
        },
        PresetDef {
            preset_id: "double_bifurcation_venturi_rect",
            preset_fn: "double_bifurcation_venturi_rect",
            family: "multilevel_composite",
            support: "parametric+ga",
            mapped_to: "D2",
            build: preset_double_bifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "triple_bifurcation_venturi_rect",
            preset_fn: "triple_bifurcation_venturi_rect",
            family: "multilevel_composite",
            support: "parametric+ga",
            mapped_to: "D3",
            build: preset_triple_bifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "double_trifurcation_venturi_rect",
            preset_fn: "double_trifurcation_venturi_rect",
            family: "mixed_tree_composite",
            support: "parametric+ga",
            mapped_to: "T9",
            build: preset_double_trifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "triple_trifurcation_venturi_rect",
            preset_fn: "triple_trifurcation_venturi_rect",
            family: "multilevel_composite",
            support: "parametric+ga",
            mapped_to: "T27",
            build: preset_triple_trifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "quad_trifurcation_venturi_rect",
            preset_fn: "quad_trifurcation_venturi_rect",
            family: "multilevel_composite",
            support: "parametric+ga",
            mapped_to: "T81",
            build: preset_quad_trifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "bifurcation_trifurcation_venturi_rect",
            preset_fn: "bifurcation_trifurcation_venturi_rect",
            family: "mixed_tree_composite",
            support: "parametric+ga",
            mapped_to: "B6",
            build: preset_bifurcation_trifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "trifurcation_bifurcation_venturi_rect",
            preset_fn: "trifurcation_bifurcation_venturi_rect",
            family: "mixed_tree_composite",
            support: "parametric+ga",
            mapped_to: "TB",
            build: preset_trifurcation_bifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "trifurcation_bifurcation_bifurcation_venturi_rect",
            preset_fn: "trifurcation_bifurcation_bifurcation_venturi_rect",
            family: "mixed_tree_composite",
            support: "parametric+ga",
            mapped_to: "TBB",
            build: preset_trifurcation_bifurcation_bifurcation_venturi_rect,
        },
        PresetDef {
            preset_id: "cell_separation_rect",
            preset_fn: "cell_separation_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "CS/WC",
            build: preset_cell_separation_rect,
        },
        PresetDef {
            preset_id: "asymmetric_bifurcation_serpentine_rect",
            preset_fn: "asymmetric_bifurcation_serpentine_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "AB",
            build: preset_asymmetric_bifurcation_serpentine_rect,
        },
        PresetDef {
            preset_id: "constriction_expansion_array_rect",
            preset_fn: "constriction_expansion_array_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "CE",
            build: preset_constriction_expansion_array_rect,
        },
        PresetDef {
            preset_id: "spiral_channel_rect",
            preset_fn: "spiral_channel_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "SP",
            build: preset_spiral_channel_rect,
        },
        PresetDef {
            preset_id: "parallel_microchannel_array_rect",
            preset_fn: "parallel_microchannel_array_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "PM",
            build: preset_parallel_microchannel_array_rect,
        },
        PresetDef {
            preset_id: "cascade_center_trifurcation_rect",
            preset_fn: "cascade_center_trifurcation_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "CCT",
            build: preset_cascade_center_trifurcation_rect,
        },
        PresetDef {
            preset_id: "incremental_filtration_tri_bi_rect",
            preset_fn: "incremental_filtration_tri_bi_rect",
            family: "specialized_composite",
            support: "parametric+ga",
            mapped_to: "CIF",
            build: preset_incremental_filtration_tri_bi_rect,
        },
    ]
}

fn inspect_preset(def: &PresetDef) -> PresetRow {
    let bp = (def.build)(def.preset_id);
    let venturi_channels = bp
        .channels
        .iter()
        .filter(|c| c.id.as_str().contains("throat"))
        .count();
    let serpentine_like_channels = bp
        .channels
        .iter()
        .filter(|c| {
            let id = c.id.as_str();
            id.contains("segment")
                || id.contains("serp")
                || id.contains("spiral")
                || id.starts_with("wide_")
                || id.starts_with("narrow_")
                || id.starts_with("ch_")
        })
        .count();

    PresetRow {
        preset_id: def.preset_id.to_string(),
        preset_fn: def.preset_fn.to_string(),
        family: def.family.to_string(),
        support: def.support.to_string(),
        mapped_to: def.mapped_to.to_string(),
        nodes: bp.nodes.len(),
        inlets: bp.inlet_count(),
        outlets: bp.outlet_count(),
        junctions: bp.junction_count(),
        channels: bp.channels.len(),
        total_length_mm: bp.total_length_m() * 1e3,
        venturi_channels,
        serpentine_like_channels,
    }
}

fn is_hydrodynamic_cavitation_ready(metrics: &SdtMetrics) -> bool {
    metrics.pressure_feasible
        && metrics.fda_main_compliant
        && metrics.cavitation_number.is_finite()
        && metrics.cavitation_number < 1.0
}

fn quantile(sorted_values: &[f64], q: f64) -> f64 {
    if sorted_values.is_empty() {
        return 0.0;
    }
    let q = q.clamp(0.0, 1.0);
    let n = sorted_values.len();
    let pos = q * (n.saturating_sub(1)) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted_values[lo]
    } else {
        let weight = pos - lo as f64;
        sorted_values[lo] * (1.0 - weight) + sorted_values[hi] * weight
    }
}

fn write_preset_catalog_csv(rows: &[PresetRow], out_path: &Path) -> Result<(), std::io::Error> {
    let mut csv = String::new();
    csv.push_str("preset_id,preset_fn,family,support,mapped_to,nodes,inlets,outlets,junctions,channels,total_length_mm,venturi_channels,serpentine_like_channels\n");
    for row in rows {
        csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{},{},{:.3},{},{}\n",
            csv_escape(&row.preset_id),
            csv_escape(&row.preset_fn),
            csv_escape(&row.family),
            csv_escape(&row.support),
            csv_escape(&row.mapped_to),
            row.nodes,
            row.inlets,
            row.outlets,
            row.junctions,
            row.channels,
            row.total_length_mm,
            row.venturi_channels,
            row.serpentine_like_channels,
        ));
    }
    std::fs::write(out_path, csv)
}

fn write_topology_mode_csv(
    rows: &[TopologyModeRow],
    out_path: &Path,
) -> Result<(), std::io::Error> {
    let mut csv = String::new();
    csv.push_str("mode_key,mode_label,topology_short,topology_name,total_candidates,evaluated_candidates,failed_metrics,feasible_count,feasible_rate_pct,hydro_cav_ready_count,best_score,median_feasible_score,p90_feasible_score,mean_feasible_score,best_candidate_id,best_sigma,best_hi_per_pass,best_bulk_hi_per_pass,best_sep3_eff,best_wbc_recovery,best_local_hct_venturi,best_cancer_dose,best_venturi_flow_frac,best_rbc_venturi_exposure,best_delta_p_kpa\n");
    for row in rows {
        csv.push_str(&format!(
            "{},{},{},{},{},{},{},{},{:.3},{},{:.6},{:.6},{:.6},{:.6},{},{:.6},{:.6e},{:.6e},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.3}\n",
            csv_escape(&row.mode_key),
            csv_escape(&row.mode_label),
            csv_escape(&row.topology_short),
            csv_escape(&row.topology_name),
            row.total_candidates,
            row.evaluated_candidates,
            row.failed_metrics,
            row.feasible_count,
            row.feasible_rate_pct,
            row.hydro_cav_ready_count,
            row.best_score,
            row.median_feasible_score,
            row.p90_feasible_score,
            row.mean_feasible_score,
            csv_escape(&row.best_candidate_id),
            row.best_sigma,
            row.best_hi_per_pass,
            row.best_bulk_hi_per_pass,
            row.best_sep3_eff,
            row.best_wbc_recovery,
            row.best_local_hct_venturi,
            row.best_cancer_dose,
            row.best_venturi_flow_frac,
            row.best_rbc_venturi_exposure,
            row.best_delta_p_kpa,
        ));
    }
    std::fs::write(out_path, csv)
}

fn write_top_candidates_csv(
    rows: &[TopCandidateRow],
    out_path: &Path,
) -> Result<(), std::io::Error> {
    let mut csv = String::new();
    csv.push_str("mode_key,mode_label,rank,candidate_id,topology_short,topology_name,score,sigma,hi_per_pass,bulk_hi_per_pass,sep3_eff,wbc_recovery,coverage_pct,delta_p_kpa,local_hct_venturi,cancer_dose_fraction,venturi_flow_fraction,rbc_venturi_exposure_fraction,fda_main,pressure_feasible,hydro_cav_ready\n");
    for row in rows {
        csv.push_str(&format!(
            "{},{},{},{},{},{},{:.6},{:.6},{:.6e},{:.6e},{:.6},{:.6},{:.3},{:.3},{:.6},{:.6},{:.6},{:.6},{},{},{}\n",
            csv_escape(&row.mode_key),
            csv_escape(&row.mode_label),
            row.rank,
            csv_escape(&row.candidate_id),
            csv_escape(&row.topology_short),
            csv_escape(&row.topology_name),
            row.score,
            row.sigma,
            row.hi_per_pass,
            row.bulk_hi_per_pass,
            row.sep3_eff,
            row.wbc_recovery,
            row.coverage_pct,
            row.delta_p_kpa,
            row.local_hct_venturi,
            row.cancer_dose_fraction,
            row.venturi_flow_fraction,
            row.rbc_venturi_exposure_fraction,
            row.fda_main,
            row.pressure_feasible,
            row.hydro_cav_ready,
        ));
    }
    std::fs::write(out_path, csv)
}

fn write_markdown_report(
    out_path: &Path,
    preset_rows: &[PresetRow],
    mode_summaries: &[ModeSummary],
    topology_rows: &[TopologyModeRow],
    top_candidates_by_mode: &BTreeMap<String, Vec<TopCandidateRow>>,
) -> Result<(), std::io::Error> {
    let mut md = String::new();
    md.push_str("# SDT Scenario Audit: cfd-schematics × cfd-optim\n\n");
    md.push_str("This report expands beyond top-5 snapshots and evaluates the full scenario surface currently available in the workspace.\n\n");
    md.push_str("## Step-by-Step Workflow\n\n");
    md.push_str("1. Enumerate every exported preset scenario from `cfd-schematics`.\n");
    md.push_str("2. Classify each preset as optimizer-mapped or schematics-only.\n");
    md.push_str("3. Run exhaustive simulation over the full `build_candidate_space()` set.\n");
    md.push_str("4. Score every evaluated candidate across eight operating modes.\n");
    md.push_str(
        "5. Rank by mode, then compare topology-level feasibility and performance spread.\n\n",
    );

    md.push_str("## 1) cfd-schematics Scenario Surface\n\n");
    md.push_str(
        "| Preset | Family | Support | Mapped To | Nodes | Channels | Outlets | Length mm |\n",
    );
    md.push_str("|---|---|---|---|---:|---:|---:|---:|\n");
    for row in preset_rows {
        md.push_str(&format!(
            "| `{}` | {} | {} | {} | {} | {} | {} | {:.2} |\n",
            row.preset_fn,
            row.family,
            row.support,
            if row.mapped_to.is_empty() {
                "—"
            } else {
                row.mapped_to.as_str()
            },
            row.nodes,
            row.channels,
            row.outlets,
            row.total_length_mm,
        ));
    }
    md.push('\n');

    let unsupported: Vec<&PresetRow> = preset_rows
        .iter()
        .filter(|r| r.support == "schematics-only")
        .collect();
    md.push_str(&format!(
        "Schematics-only presets (not in current optimizer sweep): **{}**.\n\n",
        unsupported.len()
    ));
    if !unsupported.is_empty() {
        md.push_str("Unmapped preset functions:\n");
        for row in unsupported {
            md.push_str(&format!("- `{}`\n", row.preset_fn));
        }
        md.push('\n');
    }

    md.push_str("## 2) Mode-Level Feasibility Summary\n\n");
    md.push_str(
        "| Mode | Feasible Candidates | Feasible Rate % | Best Score | Mean Feasible Score |\n",
    );
    md.push_str("|---|---:|---:|---:|---:|\n");
    for summary in mode_summaries {
        md.push_str(&format!(
            "| {} | {} | {:.2} | {:.4} | {:.4} |\n",
            summary.mode_label,
            summary.feasible_candidates,
            summary.feasible_rate_pct,
            summary.best_score,
            summary.mean_feasible_score,
        ));
    }
    md.push('\n');

    md.push_str("## 3) Detailed Simulation Results by Mode\n\n");
    for summary in mode_summaries {
        md.push_str(&format!("### {}\n\n", summary.mode_label));

        let mut mode_topology_rows: Vec<&TopologyModeRow> = topology_rows
            .iter()
            .filter(|r| r.mode_key == summary.mode_key)
            .collect();
        mode_topology_rows.sort_by(|a, b| b.best_score.total_cmp(&a.best_score));

        md.push_str("Top topology families (by best feasible score):\n\n");
        md.push_str(
            "| Topology | Total | Feasible | Rate % | Best | Median | P90 | Hydro+Cav Ready |\n",
        );
        md.push_str("|---|---:|---:|---:|---:|---:|---:|---:|\n");
        for row in mode_topology_rows.iter().take(10) {
            md.push_str(&format!(
                "| {} ({}) | {} | {} | {:.1} | {:.4} | {:.4} | {:.4} | {} |\n",
                row.topology_name,
                row.topology_short,
                row.total_candidates,
                row.feasible_count,
                row.feasible_rate_pct,
                row.best_score,
                row.median_feasible_score,
                row.p90_feasible_score,
                row.hydro_cav_ready_count,
            ));
        }
        md.push('\n');

        if let Some(top_rows) = top_candidates_by_mode.get(&summary.mode_key) {
            md.push_str("Top candidate designs:\n\n");
            md.push_str("| Rank | Candidate ID | Topology | Score | sigma | HI/pass | Bulk HI | Local Hct@V | Cancer Dose | RBC Venturi | Venturi Qfrac | Sep3 | WBC rec | dP kPa | Hydro+Cav |\n");
            md.push_str(
                "|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|\n",
            );
            for row in top_rows.iter().take(10) {
                let sigma = if row.sigma.is_finite() {
                    format!("{:.4}", row.sigma)
                } else {
                    "inf".to_string()
                };
                md.push_str(&format!(
                    "| {} | `{}` | {} ({}) | {:.4} | {} | {:.2e} | {:.2e} | {:.3} | {:.3} | {:.3} | {:.3} | {:.4} | {:.4} | {:.2} | {} |\n",
                    row.rank,
                    row.candidate_id,
                    row.topology_name,
                    row.topology_short,
                    row.score,
                    sigma,
                    row.hi_per_pass,
                    row.bulk_hi_per_pass,
                    row.local_hct_venturi,
                    row.cancer_dose_fraction,
                    row.rbc_venturi_exposure_fraction,
                    row.venturi_flow_fraction,
                    row.sep3_eff,
                    row.wbc_recovery,
                    row.delta_p_kpa,
                    if row.hydro_cav_ready { "PASS" } else { "FAIL" },
                ));
            }
            md.push('\n');
        }
    }

    md.push_str("## 4) Interpretation\n\n");
    md.push_str("- The optimizer evaluates a broad multi-topology candidate set, but only a subset of cfd-schematics presets are presently mapped into the sweep.\n");
    md.push_str("- Hydrodynamic-cavitation-ready designs are mode-dependent; ranking changes when separation-heavy objectives are used.\n");
    md.push_str("- The GA-only `AdaptiveTree` scenario remains available through `DesignCandidate::to_blueprint`, but is not part of `build_candidate_space()`.\n");

    std::fs::write(out_path, md)
}

fn write_cif_cct_stepwise_markdown(
    out_path: &Path,
    combined_ranked: &[(f64, &EvaluatedRecord)],
) -> Result<(), std::io::Error> {
    if let Some(parent) = out_path.parent() {
        std::fs::create_dir_all(parent)?;
    }

    const MAX_ROWS_PER_TOPOLOGY: usize = 30;

    let mut total_cif = 0usize;
    let mut total_cct = 0usize;
    for (_, rec) in combined_ranked {
        match rec.candidate.topology {
            cfd_optim::DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => total_cif += 1,
            cfd_optim::DesignTopology::CascadeCenterTrifurcationSeparator { .. } => total_cct += 1,
            _ => {}
        }
    }

    let mut md = String::new();
    md.push_str("# Milestone 12 Scenario Audit Appendix\n\n");
    md.push_str("## Combined SDT + Leukapheresis: Stepwise CIF/CCT Analysis\n\n");
    md.push_str("This appendix extracts only the feasible Combined-mode candidates and reports staged skimming controls for CIF and CCT designs.\n\n");
    md.push_str(&format!(
        "- Feasible combined candidates: **{}**\n- CCT candidates in feasible set: **{}**\n- CIF candidates in feasible set: **{}**\n\n",
        combined_ranked.len(),
        total_cct,
        total_cif
    ));
    md.push_str("Skimming-factor definitions used below:\n");
    md.push_str("- `q_center(frac) = frac^3 / (frac^3 + 2*((1-frac)/2)^3)`\n");
    md.push_str(
        "- CIF model venturi-flow fraction: `q_pretri^n_pretri * q_terminal_tri * q_bi_treat`\n",
    );
    md.push_str("- CCT model venturi-flow fraction: `q_center^n_levels`\n\n");

    md.push_str("## CIF: Controlled Incremental Filtration (Trifurcation -> Trifurcation -> Bifurcation)\n\n");
    md.push_str("| Global rank | Candidate ID | Score | n_pretri | pcf | tcf | btf | q_pretri | q_terminal_tri | q_bi_treat | model_qfrac | solved_qfrac | q_pretri_solved | q_tri_solved | q_bi_solved | tail_mm | cancer_targeted_cav | wbc_recovery | rbc_venturi_exposure | HI/pass | HI15m(3kg)% | ECV mL |\n");
    md.push_str(
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n",
    );

    let mut cif_written = 0usize;
    for (idx, (score, rec)) in combined_ranked.iter().enumerate() {
        let n_pretri = match rec.candidate.topology {
            cfd_optim::DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } => n_pretri,
            _ => continue,
        };
        if cif_written >= MAX_ROWS_PER_TOPOLOGY {
            break;
        }

        let pcf = rec.candidate.cif_pretri_center_frac();
        let tcf = rec.candidate.cif_terminal_tri_center_frac();
        let btf = rec.candidate.cif_terminal_bi_treat_frac();
        let q_pretri = cfd_1d::cell_separation::tri_center_q_frac(pcf);
        let q_terminal_tri = cfd_1d::cell_separation::tri_center_q_frac(tcf);
        let model_qfrac =
            (q_pretri.powi(i32::from(n_pretri)) * q_terminal_tri * btf).clamp(0.0, 1.0);

        md.push_str(&format!(
            "| {} | `{}` | {:.4} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.2} | {:.4} | {:.4} | {:.4} | {:.2e} | {:.3} | {:.3} |\n",
            idx + 1,
            rec.candidate.id,
            score,
            n_pretri,
            pcf,
            tcf,
            btf,
            q_pretri,
            q_terminal_tri,
            btf,
            model_qfrac,
            rec.metrics.venturi_flow_fraction,
            rec.metrics.cif_pretri_qfrac_mean,
            rec.metrics.cif_terminal_tri_qfrac,
            rec.metrics.cif_terminal_bi_qfrac,
            rec.metrics.cif_outlet_tail_length_mm,
            rec.metrics.cancer_targeted_cavitation,
            rec.metrics.wbc_recovery,
            rec.metrics.rbc_venturi_exposure_fraction,
            rec.metrics.hemolysis_index_per_pass,
            rec.metrics.projected_hemolysis_15min_pediatric_3kg * 100.0,
            rec.metrics.total_ecv_ml,
        ));
        cif_written += 1;
    }
    if cif_written == 0 {
        md.push_str("| - | _No feasible CIF candidate in combined-mode set_ | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - | - |\n");
    }
    md.push('\n');

    md.push_str("## CCT: Cascade Center Trifurcation\n\n");
    md.push_str("| Global rank | Candidate ID | Score | n_levels | center_frac | q_center | model_qfrac | solved_qfrac | q_center_solved_mean | cancer_targeted_cav | wbc_recovery | rbc_venturi_exposure | HI/pass | HI15m(3kg)% | ECV mL |\n");
    md.push_str("|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    let mut cct_written = 0usize;
    for (idx, (score, rec)) in combined_ranked.iter().enumerate() {
        let n_levels = match rec.candidate.topology {
            cfd_optim::DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } => n_levels,
            _ => continue,
        };
        if cct_written >= MAX_ROWS_PER_TOPOLOGY {
            break;
        }

        let center_frac = rec.candidate.trifurcation_center_frac.clamp(0.20, 0.70);
        let q_center = cfd_1d::cell_separation::tri_center_q_frac(center_frac);
        let model_qfrac = q_center.powi(i32::from(n_levels)).clamp(0.0, 1.0);

        md.push_str(&format!(
            "| {} | `{}` | {:.4} | {} | {:.3} | {:.3} | {:.3} | {:.3} | {:.3} | {:.4} | {:.4} | {:.4} | {:.2e} | {:.3} | {:.3} |\n",
            idx + 1,
            rec.candidate.id,
            score,
            n_levels,
            center_frac,
            q_center,
            model_qfrac,
            rec.metrics.venturi_flow_fraction,
            rec.metrics.cct_stage_center_qfrac_mean,
            rec.metrics.cancer_targeted_cavitation,
            rec.metrics.wbc_recovery,
            rec.metrics.rbc_venturi_exposure_fraction,
            rec.metrics.hemolysis_index_per_pass,
            rec.metrics.projected_hemolysis_15min_pediatric_3kg * 100.0,
            rec.metrics.total_ecv_ml,
        ));
        cct_written += 1;
    }
    if cct_written == 0 {
        md.push_str("| - | _No feasible CCT candidate in combined-mode set_ | - | - | - | - | - | - | - | - | - | - | - | - | - |\n");
    }
    md.push('\n');

    md.push_str("## Observations\n\n");
    md.push_str("- CIF entries expose separate control knobs for pre-cascade center bias (`pcf`), terminal trifurcation center bias (`tcf`), and terminal bifurcation treatment routing (`btf`).\n");
    md.push_str("- CCT entries isolate the effect of repeated center-only trifurcation via (`n_levels`, `center_frac`).\n");
    md.push_str("- Divergence between `model_qfrac` and `solved_qfrac` indicates downstream resistance coupling that the local split model does not capture.\n");
    md.push_str("- `HI15m(3kg)%` reports projected cumulative hemolysis over a 15-minute pediatric (3 kg) treatment window.\n");
    md.push_str("- Candidate ranking in this appendix is the global ranking from Combined SDT + leukapheresis mode.\n");

    std::fs::write(out_path, md)
}

fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

// ── Canonical preset builders used by the surface audit ─────────────────────

fn preset_venturi_chain(name: &str) -> NetworkBlueprint {
    venturi_chain(name, 45e-3, 4e-3, 150e-6)
}

fn preset_venturi_rect(name: &str) -> NetworkBlueprint {
    venturi_rect(name, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_symmetric_bifurcation(name: &str) -> NetworkBlueprint {
    symmetric_bifurcation(name, 12e-3, 6e-3, 2e-3, 1e-3)
}

fn preset_bifurcation_rect(name: &str) -> NetworkBlueprint {
    bifurcation_rect(name, 12e-3, 6e-3, 2e-3, 1e-3, 1e-3)
}

fn preset_symmetric_trifurcation(name: &str) -> NetworkBlueprint {
    symmetric_trifurcation(name, 12e-3, 6e-3, 2e-3, 1e-3)
}

fn preset_trifurcation_rect(name: &str) -> NetworkBlueprint {
    trifurcation_rect(name, 12e-3, 6e-3, 2e-3, 1e-3, 1e-3)
}

fn preset_serpentine_chain(name: &str) -> NetworkBlueprint {
    serpentine_chain(name, 6, 7.5e-3, 2e-3)
}

fn preset_serpentine_rect(name: &str) -> NetworkBlueprint {
    serpentine_rect(name, 6, 7.5e-3, 2.0e-3, 1.0e-3)
}

fn preset_venturi_serpentine_rect(name: &str) -> NetworkBlueprint {
    venturi_serpentine_rect(name, 2.0e-3, 100e-6, 1.0e-3, 300e-6, 6, 7.5e-3)
}

fn preset_serial_double_venturi_rect(name: &str) -> NetworkBlueprint {
    serial_double_venturi_rect(name, 2.0e-3, 100e-6, 1.0e-3, 300e-6, 7.5e-3)
}

fn preset_bifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    bifurcation_venturi_rect(name, 15e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_bifurcation_serpentine_rect(name: &str) -> NetworkBlueprint {
    bifurcation_serpentine_rect(name, 15e-3, 6, 7.5e-3, 2.0e-3, 1.0e-3)
}

fn preset_trifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    trifurcation_venturi_rect(name, 15e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_trifurcation_serpentine_rect(name: &str) -> NetworkBlueprint {
    trifurcation_serpentine_rect(name, 15e-3, 6, 7.5e-3, 2.0e-3, 1.0e-3)
}

fn preset_double_bifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    double_bifurcation_venturi_rect(name, 12e-3, 9e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_triple_bifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    triple_bifurcation_venturi_rect(name, 10e-3, 8e-3, 6e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_double_trifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    double_trifurcation_venturi_rect(name, 12e-3, 8e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_triple_trifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    triple_trifurcation_venturi_rect(name, 9e-3, 7e-3, 5e-3, 2.0e-3, 0.45, 100e-6, 1.0e-3, 300e-6)
}

fn preset_quad_trifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    quad_trifurcation_venturi_rect(
        name, 8e-3, 6e-3, 5e-3, 4e-3, 2.0e-3, 0.45, 100e-6, 1.0e-3, 300e-6,
    )
}

fn preset_bifurcation_trifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    bifurcation_trifurcation_venturi_rect(name, 12e-3, 8e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_trifurcation_bifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    trifurcation_bifurcation_venturi_rect(name, 12e-3, 8e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6)
}

fn preset_trifurcation_bifurcation_bifurcation_venturi_rect(name: &str) -> NetworkBlueprint {
    trifurcation_bifurcation_bifurcation_venturi_rect(
        name, 10e-3, 7e-3, 5e-3, 2.0e-3, 0.45, 100e-6, 1.0e-3, 300e-6,
    )
}

fn preset_cell_separation_rect(name: &str) -> NetworkBlueprint {
    cell_separation_rect(name, 22.5e-3, 2.0e-3, 100e-6, 1.0e-3, 300e-6, 22.5e-3)
}

fn preset_asymmetric_bifurcation_serpentine_rect(name: &str) -> NetworkBlueprint {
    asymmetric_bifurcation_serpentine_rect(name, 15e-3, 6, 7.5e-3, 0.5, 2.0e-3, 1.0e-3)
}

fn preset_constriction_expansion_array_rect(name: &str) -> NetworkBlueprint {
    constriction_expansion_array_rect(name, 10, 3.0e-3, 1.5e-3, 200e-6, 80e-6, 60e-6)
}

fn preset_spiral_channel_rect(name: &str) -> NetworkBlueprint {
    spiral_channel_rect(name, 8, 5.0e-3, 200e-6, 60e-6)
}

fn preset_parallel_microchannel_array_rect(name: &str) -> NetworkBlueprint {
    parallel_microchannel_array_rect(name, 64, 45e-3, 200e-6, 60e-6)
}

fn preset_cascade_center_trifurcation_rect(name: &str) -> NetworkBlueprint {
    cascade_center_trifurcation_rect(name, 12e-3, 8e-3, 2, 2.0e-3, 0.45, 100e-6, 300e-6, 1.0e-3)
}

fn preset_incremental_filtration_tri_bi_rect(name: &str) -> NetworkBlueprint {
    incremental_filtration_tri_bi_rect(
        name, 12e-3, 8e-3, 6e-3, 2, 2.0e-3, 0.45, 0.68, 100e-6, 300e-6, 1.0e-3,
    )
}
