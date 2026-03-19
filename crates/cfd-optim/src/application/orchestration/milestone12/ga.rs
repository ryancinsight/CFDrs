use crate::application::orchestration::{
    ensure_release_reports, fast_env, fast_mode, ga_matches_lineage_sequence, init_tracing,
    is_selective_report_topology, milestone12_ranked_pool_size, resolve_output_directories,
    save_figure,
};
use crate::application::search::mutations::{
    promote_option1_candidate_to_ga_seed, promote_option2_candidate_to_ga_seed,
};
use crate::application::search::genetic::BlueprintGeneticOptimizer;
use crate::delivery::{load_top5_report_json, save_pareto_points, save_top5_report_json};
use crate::domain::OptimizationGoal;
use crate::reporting::{
    is_hydrosdt_venturi_report_candidate, pareto_pool_from_report_designs,
    rank_ga_hydrosdt_report_designs,
    milestone12_lineage_key, validate_milestone12_candidate, write_goal_audit_report,
    GoalAuditEntry, GoalAuditStatus, Milestone12LineageKey, Milestone12ReportDesign,
    Milestone12Stage, ParetoPoint, ParetoTag,
};


use super::report::{write_stage_summary, Milestone12GaSummary, GA_SUMMARY_PATH};
use super::types::{Milestone12GaRun, Milestone12StageArtifact};

#[derive(Debug, Clone, serde::Serialize)]
struct GaLineageFamilyCount {
    family: String,
    lane: String,
    count: usize,
}

#[derive(Debug, Clone, serde::Serialize)]
struct GaLineageAuditEntry {
    rank: usize,
    candidate_id: String,
    score: f64,
    adjusted_selection_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    dominant_family: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    dominant_lane: Option<String>,
    geometry_concentration_penalty: f64,
    operating_point_diversity_penalty: f64,
    serpentine_repeat_penalty: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    family_counts: Vec<GaLineageFamilyCount>,
}

#[derive(Debug, Clone, serde::Serialize)]
struct GaLineageAncestryEntry {
    rank: usize,
    candidate_id: String,
    score: f64,
    lineage_depth: usize,
    total_inherited_mutation_count: usize,
    geometry_event_count: usize,
    operating_point_event_count: usize,
}

fn parse_lineage_family_lane(metadata: &str) -> Option<(String, String)> {
    let mut family = None;
    let mut lane = None;
    for part in metadata.split(';') {
        let (key, value) = part.split_once('=')?;
        match key {
            "family" => family = Some(value.to_string()),
            "lane" => lane = Some(value.to_string()),
            _ => {}
        }
    }
    Some((family?, lane?))
}

fn parse_lineage_operator(metadata: &str) -> Option<String> {
    metadata.split(';').find_map(|part| {
        let (key, value) = part.split_once('=')?;
        (key == "operator").then(|| value.to_string())
    })
}

fn positive_ratio(value: f64) -> f64 {
    value.max(0.0)
}

fn lane_serpentine_diminishing_return(
    route: &cfd_schematics::ChannelRouteSpec,
) -> f64 {
    let Some(serpentine) = route.serpentine.as_ref() else {
        return 0.0;
    };

    let neutral_segments = 4.0;
    let neutral_bend_radius_m = (route.width_m * 2.5).max(route.height_m);
    let neutral_segment_length_m = (route.length_m / 4.0).max(route.width_m);

    let segments_factor = serpentine.segments as f64 / neutral_segments;
    let bend_factor = serpentine.bend_radius_m / neutral_bend_radius_m.max(1.0e-12);
    let length_factor = serpentine.segment_length_m / neutral_segment_length_m.max(1.0e-12);

    let compact_depth = positive_ratio((1.0 - bend_factor) / 0.20)
        + positive_ratio((segments_factor - 1.25) / 0.50)
        + positive_ratio((1.0 - length_factor) / 0.20);
    let dense_depth = positive_ratio((segments_factor - 1.75) / 0.75)
        + positive_ratio((0.85 - length_factor) / 0.15)
        + positive_ratio((0.85 - bend_factor) / 0.15);
    let smooth_depth = positive_ratio((bend_factor - 1.30) / 0.30)
        + positive_ratio((length_factor - 1.10) / 0.25);
    let long_depth = positive_ratio((length_factor - 1.30) / 0.30)
        + positive_ratio((segments_factor - 1.25) / 0.75);

    let dominant_family_depth = compact_depth
        .max(dense_depth)
        .max(smooth_depth)
        .max(long_depth)
        .max(0.0);
    (dominant_family_depth - 1.0).max(0.0)
}

fn serpentine_repeat_penalty(design: &Milestone12ReportDesign) -> f64 {
    let topology = match design.candidate.topology_spec() {
        Ok(topology) => topology,
        Err(_) => return 0.0,
    };
    let mut lane_family_counts: std::collections::HashMap<(String, String), usize> =
        std::collections::HashMap::new();
    let mut geometry_event_count = 0_usize;
    if let Some(lineage) = design.candidate.blueprint().lineage() {
        for event in &lineage.mutations {
            let operator = parse_lineage_operator(&event.mutation).unwrap_or_default();
            if !operator.starts_with("operating_point") {
                geometry_event_count += 1;
            }
            if let Some((family, lane)) = parse_lineage_family_lane(&event.mutation) {
                if family.starts_with("serpentine_") {
                    *lane_family_counts.entry((lane, family)).or_insert(0) += 1;
                }
            }
        }
    }

    let geometry_depth_scale = (geometry_event_count as f64).ln_1p().max(1.0);

    lane_family_counts
        .into_iter()
        .filter_map(|((lane, family), count)| {
            topology.channel_route(&lane).map(|route| {
                let repeat_depth = count.saturating_sub(1) as f64;
                let ancestry_share = count as f64 / geometry_event_count.max(1) as f64;
                let family_bias = if family == "serpentine_compact" { 1.2 } else { 1.0 };
                repeat_depth
                    * lane_serpentine_diminishing_return(route)
                    * family_bias
                    * (1.0 + ancestry_share * geometry_depth_scale)
            })
        })
        .fold(0.0_f64, f64::max)
}

fn operating_point_diversity_penalty(design: &Milestone12ReportDesign) -> f64 {
    let mut operating_point_family_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    if let Some(lineage) = design.candidate.blueprint().lineage() {
        for event in &lineage.mutations {
            let operator = parse_lineage_operator(&event.mutation).unwrap_or_default();
            if !operator.starts_with("operating_point") {
                continue;
            }
            if let Some((family, _lane)) = parse_lineage_family_lane(&event.mutation) {
                *operating_point_family_counts.entry(family).or_insert(0) += 1;
            }
        }
    }

    let total_events: usize = operating_point_family_counts.values().sum();
    if total_events <= 1 {
        return 0.0;
    }

    let unique_families = operating_point_family_counts.len();
    let dominant_count = operating_point_family_counts
        .values()
        .copied()
        .max()
        .unwrap_or(0);
    let repeated_events = total_events.saturating_sub(unique_families) as f64;
    let dominance_share = dominant_count as f64 / total_events as f64;
    repeated_events * (1.0 + dominance_share)
}

fn build_ga_lineage_audit_entries(
    ranked: &[Milestone12ReportDesign],
    candidate_scores: &std::collections::HashMap<(String, u64), f64>,
) -> Vec<GaLineageAuditEntry> {
    ranked
        .iter()
        .map(|design| {
            let mut counts: std::collections::HashMap<(String, String), usize> =
                std::collections::HashMap::new();
            if let Some(lineage) = design.candidate.blueprint().lineage() {
                for event in &lineage.mutations {
                    if let Some((family, lane)) = parse_lineage_family_lane(&event.mutation) {
                        *counts.entry((family, lane)).or_insert(0) += 1;
                    }
                }
            }
            let mut family_counts: Vec<GaLineageFamilyCount> = counts
                .into_iter()
                .map(|((family, lane), count)| GaLineageFamilyCount { family, lane, count })
                .collect();
            family_counts.sort_by(|left, right| {
                right
                    .count
                    .cmp(&left.count)
                    .then_with(|| left.lane.cmp(&right.lane))
                    .then_with(|| left.family.cmp(&right.family))
            });
            let dominant_family = family_counts.first().map(|entry| entry.family.clone());
            let dominant_lane = family_counts.first().map(|entry| entry.lane.clone());
            let geometry_penalty = serpentine_repeat_penalty(design);
            let operating_point_penalty = operating_point_diversity_penalty(design);
            GaLineageAuditEntry {
                rank: design.rank,
                candidate_id: design.candidate.id.clone(),
                score: design.score,
                adjusted_selection_score: candidate_scores
                    .get(&(design.candidate.id.clone(), design.score.to_bits()))
                    .copied()
                    .unwrap_or(design.score),
                dominant_family,
                dominant_lane,
                geometry_concentration_penalty: geometry_penalty,
                operating_point_diversity_penalty: operating_point_penalty,
                serpentine_repeat_penalty: geometry_penalty,
                family_counts,
            }
        })
        .collect()
}

fn build_ga_lineage_ancestry_entries(
    ranked: &[Milestone12ReportDesign],
) -> Vec<GaLineageAncestryEntry> {
    ranked
        .iter()
        .map(|design| {
            let lineage_depth = design
                .candidate
                .blueprint()
                .lineage()
                .map_or(0, |lineage| lineage.mutations.len());
            let (geometry_event_count, operating_point_event_count) = design
                .candidate
                .blueprint()
                .lineage()
                .map(|lineage| {
                    lineage.mutations.iter().fold((0_usize, 0_usize), |mut counts, event| {
                        let operator = parse_lineage_operator(&event.mutation).unwrap_or_default();
                        if operator.starts_with("operating_point") {
                            counts.1 += 1;
                        } else {
                            counts.0 += 1;
                        }
                        counts
                    })
                })
                .unwrap_or((0, 0));
            GaLineageAncestryEntry {
                rank: design.rank,
                candidate_id: design.candidate.id.clone(),
                score: design.score,
                lineage_depth,
                total_inherited_mutation_count: lineage_depth.saturating_sub(1),
                geometry_event_count,
                operating_point_event_count,
            }
        })
        .collect()
}

pub fn run_milestone12_ga() -> Result<Milestone12GaRun, Box<dyn std::error::Error>> {
    init_tracing();
    ensure_release_reports()?;
    let (_, out_dir, figures_dir) = resolve_output_directories()?;
    let is_fast = fast_mode();
    let retained_pool_size = milestone12_ranked_pool_size();
    let ga_population = if is_fast {
        fast_env("M12_FAST_GA_POPULATION", 50)
    } else {
        80
    };
    let ga_generations = if is_fast {
        fast_env("M12_FAST_GA_GENERATIONS", 16)
    } else {
        120
    };

    // Load seeds from previous-stage disk artifacts instead of rebuilding
    // the entire 80K candidate space (~1.2 GB). The pipeline guarantees
    // Option 1 and Option 2 run before GA.
    let option1_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option1_ultrasound_top5.json"))?;
    let option2_ranked =
        load_top5_report_json(&out_dir.join("two_concept_option2_venturi_top5.json"))?;

    let mut seeds = Vec::with_capacity(6);
    for design in option1_ranked.iter().take(3) {
        if let Ok(promoted) = promote_option1_candidate_to_ga_seed(&design.candidate) {
            seeds.push(promoted);
        }
    }
    for design in option2_ranked.iter().take(3) {
        if let Ok(promoted) = promote_option2_candidate_to_ga_seed(&design.candidate) {
            seeds.push(promoted);
        }
    }
    if seeds.is_empty() {
        return Err("GA seed selection: no Option 1 or Option 2 ranked designs on disk".into());
    }
    let option2_baseline = option2_ranked
        .first()
        .cloned()
        .ok_or("GA seed selection: no Option 2 ranked designs on disk")?;
    drop(option1_ranked); // free the full designs
    let selected_lineage_key: Option<Milestone12LineageKey> = option2_ranked
        .first()
        .and_then(|design| milestone12_lineage_key(&design.candidate));

    let ga_result =
        BlueprintGeneticOptimizer::new(OptimizationGoal::InPlaceDeanSerpentineRefinement)
            .with_seeds(seeds)
            .with_hydrosdt_baseline_serial_cavitation_dose_fraction(
                option2_baseline.metrics.serial_cavitation_dose_fraction,
            )
            .with_population(ga_population)
            .with_max_generations(ga_generations)
            .with_top_k(retained_pool_size)
            .run()?;

    let mut candidate_order: Vec<_> = ga_result
        .all_candidates
        .iter()
        .filter(|ranked_candidate| is_selective_report_topology(&ranked_candidate.candidate))
        .collect();
    candidate_order.sort_unstable_by(|left, right| {
        let left_exceeds = left.evaluation.exceeds_all_baselines.unwrap_or(false);
        let right_exceeds = right.evaluation.exceeds_all_baselines.unwrap_or(false);
        let left_lineage = selected_lineage_key
            .as_ref()
            .is_some_and(|lineage_key| ga_matches_lineage_sequence(&left.candidate, lineage_key));
        let right_lineage = selected_lineage_key
            .as_ref()
            .is_some_and(|lineage_key| ga_matches_lineage_sequence(&right.candidate, lineage_key));
        right_exceeds
            .cmp(&left_exceeds)
            .then_with(|| right_lineage.cmp(&left_lineage))
            .then_with(|| right.selection_score.total_cmp(&left.selection_score))
            .then_with(|| right.evaluation.score_or_zero().total_cmp(&left.evaluation.score_or_zero()))
            .then_with(|| left.candidate.id.cmp(&right.candidate.id))
    });

    let mut ranked: Vec<Milestone12ReportDesign> = candidate_order
        .into_iter()
        .filter_map(|ranked_candidate| {
            let score = ranked_candidate.evaluation.score?;
            let design = Milestone12ReportDesign::from_blueprint_candidate(
                0,
                ranked_candidate.candidate.clone(),
                score,
            )
            .ok()?;
            validate_milestone12_candidate(&design.candidate, Milestone12Stage::GaRefined)
                .ok()?;
            Some(design)
        })
        .take(retained_pool_size)
        .collect();

    for (index, design) in ranked.iter_mut().enumerate() {
        design.rank = index + 1;
    }
    let full_hydrosdt_ranked: Vec<Milestone12ReportDesign> = ranked
        .iter()
        .filter(|design| is_hydrosdt_venturi_report_candidate(design))
        .cloned()
        .collect();
    let ranked = rank_ga_hydrosdt_report_designs(&ranked, &option2_baseline);
    if ranked.is_empty() {
        return Err("GA produced no viable designs after validation".into());
    }

    let audit_entries: Vec<GoalAuditEntry> = ga_result
        .all_candidates
        .iter()
        .map(|ranked_candidate| GoalAuditEntry {
            candidate_id: ranked_candidate.candidate.id.clone(),
            blueprint_name: ranked_candidate.candidate.blueprint.name.clone(),
            goal: OptimizationGoal::InPlaceDeanSerpentineRefinement,
            status: match ranked_candidate.evaluation.status {
                crate::application::objectives::BlueprintEvaluationStatus::Eligible => {
                    GoalAuditStatus::Eligible
                }
                crate::application::objectives::BlueprintEvaluationStatus::ScreenedOut => {
                    GoalAuditStatus::ScreenedOut
                }
            },
            reasons: ranked_candidate.evaluation.screening_reasons.clone(),
            score: ranked_candidate.evaluation.score,
            baseline_scores: ranked_candidate.evaluation.baseline_scores.clone(),
            exceeds_all_baselines: ranked_candidate.evaluation.exceeds_all_baselines,
            treatment_residence_time_s: Some(
                ranked_candidate
                    .evaluation
                    .residence
                    .treatment_residence_time_s,
            ),
            separation_efficiency: Some(
                ranked_candidate.evaluation.separation.separation_efficiency,
            ),
            cavitation_selectivity_score: Some(
                ranked_candidate.evaluation.venturi.cavitation_selectivity_score,
            ),
            peak_dean_number: Some(
                ranked_candidate
                    .evaluation
                    .venturi
                    .placements
                    .iter()
                    .map(|placement| placement.dean_number)
                    .fold(0.0_f64, f64::max),
            ),
        })
        .collect();
    let audit = write_goal_audit_report(&out_dir, "ga_audit", &audit_entries)?;
    let candidate_selection_scores: std::collections::HashMap<(String, u64), f64> = ga_result
        .all_candidates
        .iter()
        .map(|candidate| {
            (
                (
                    candidate.candidate.id.clone(),
                    candidate.evaluation.score_or_zero().to_bits(),
                ),
                candidate.selection_score,
            )
        })
        .collect();

    // Persist lightweight Pareto points for GA pool (~32 bytes each).
    // Extract TTCI/RBC-protection directly from evaluation data — avoids
    // cloning 200 full BlueprintCandidates and computing SdtMetrics for each.
    let ga_pareto_points: Vec<ParetoPoint> =
        pareto_pool_from_report_designs(&full_hydrosdt_ranked, ParetoTag::Ga, 200);
    let ga_pool_all_path = out_dir.join("ga_pool_all.json");
    save_pareto_points(&ga_pareto_points, &ga_pool_all_path)?;
    drop(ga_pareto_points);

    let top5_path = out_dir.join("ga_hydrosdt_top5.json");
    save_top5_report_json(&ranked, &top5_path)?;
    let lineage_audit_path = out_dir.join("ga_lineage_audit_top5.json");
    std::fs::write(
        &lineage_audit_path,
        serde_json::to_string_pretty(&build_ga_lineage_audit_entries(&ranked, &candidate_selection_scores))?,
    )?;
    let ancestry_audit_path = out_dir.join("ga_lineage_ancestry_top5.json");
    std::fs::write(
        &ancestry_audit_path,
        serde_json::to_string_pretty(&build_ga_lineage_ancestry_entries(&ranked))?,
    )?;
    let figure_path = figures_dir.join("top_hydrosdt_schematic.svg");
    save_figure(
        ranked[0].candidate.blueprint(),
        &figure_path,
        "Figure 6 (HydroSDT GA rank-1)",
    );
    write_stage_summary(
        &out_dir,
        GA_SUMMARY_PATH,
        &Milestone12GaSummary {
            evaluated_count: ga_result.all_candidates.len(),
            authoritative_run: true,
            fast_mode: crate::application::orchestration::fast_mode(),
            best_per_generation: ga_result.best_per_generation.clone(),
        },
    )?;

    Ok(Milestone12GaRun {
        ranked,
        best_per_generation: ga_result.best_per_generation,
        audit,
        artifacts: vec![
            Milestone12StageArtifact {
                label: "ga_top5".to_string(),
                path: top5_path,
            },
            Milestone12StageArtifact {
                label: "ga_figure".to_string(),
                path: figure_path,
            },
            Milestone12StageArtifact {
                label: "ga_lineage_audit_top5".to_string(),
                path: lineage_audit_path,
            },
            Milestone12StageArtifact {
                label: "ga_lineage_ancestry_top5".to_string(),
                path: ancestry_audit_path,
            },
        ],
    })
}
