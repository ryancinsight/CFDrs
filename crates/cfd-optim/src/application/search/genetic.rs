use std::collections::{hash_map::Entry, HashMap, HashSet};
use std::hash::{Hash, Hasher};
use std::sync::OnceLock;

use crate::application::objectives::{evaluate_goal, BlueprintObjectiveEvaluation};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::healthy_cell_protection_index;
use cfd_schematics::{TopologyLineageEvent, TopologyLineageMetadata};
use serde::{Deserialize, Serialize};

pub use super::mutations::{generate_ga_crossover_children, generate_ga_mutations};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintRankedCandidate {
    pub rank: usize,
    pub candidate: BlueprintCandidate,
    pub evaluation: BlueprintObjectiveEvaluation,
    pub selection_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintEvolutionResult {
    pub top_candidates: Vec<BlueprintRankedCandidate>,
    pub all_candidates: Vec<BlueprintRankedCandidate>,
    pub best_per_generation: Vec<f64>,
}

pub struct BlueprintGeneticOptimizer {
    goal: OptimizationGoal,
    population: usize,
    max_generations: usize,
    top_k: usize,
    seeds: Vec<BlueprintCandidate>,
    hydrosdt_baseline_serial_cavitation_dose_fraction: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
struct HydroShortlistPolicy {
    baseline_serial_cavitation_dose_fraction: f64,
    required_gain_margin: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CandidateFingerprint {
    topology_hash: u64,
    flow_rate_quantized: i64,
    inlet_gauge_quantized: i64,
    hematocrit_quantized: i64,
}

#[derive(Debug, Clone, Copy)]
struct DiversitySignature {
    flow_rate_m3_s: f64,
    inlet_gauge_pa: f64,
    throat_width_m: f64,
    serial_venturi_stages: f64,
    parallel_venturi_count: f64,
    split_stage_count: f64,
    serpentine_arm_count: f64,
    cavitation_selectivity: f64,
    separation_efficiency: f64,
    healthy_cell_protection: f64,
    residence_time_s: f64,
    dean_norm: f64,
}

fn evaluation_score_or_zero(evaluation: &BlueprintObjectiveEvaluation) -> f64 {
    evaluation.score.unwrap_or(0.0)
}

fn geometry_penalty_weight() -> f64 {
    static WEIGHT: OnceLock<f64> = OnceLock::new();
    *WEIGHT.get_or_init(|| {
        std::env::var("M12_GA_GEOMETRY_PENALTY_WEIGHT")
            .ok()
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| value.is_finite() && *value >= 0.0)
            .unwrap_or(0.00043)
    })
}

fn operating_point_penalty_weight() -> f64 {
    0.00020
}

fn hydrosdt_cavitation_gain_margin() -> f64 {
    static MARGIN: OnceLock<f64> = OnceLock::new();
    *MARGIN.get_or_init(|| {
        std::env::var("M12_GA_CAVITATION_GAIN_MARGIN")
            .ok()
            .and_then(|value| value.parse::<f64>().ok())
            .filter(|value| value.is_finite() && *value >= 0.0)
            .unwrap_or(0.01)
    })
}

fn cavitation_strength_from_sigma(cavitation_number: f64) -> f64 {
    let strength = (1.0 - cavitation_number).max(0.0);
    strength / (1.0 + strength)
}

fn hydrosdt_serial_cavitation_dose_fraction(candidate: &BlueprintRankedCandidate) -> Option<f64> {
    let topology = candidate.candidate.topology_spec().ok()?;
    let serial_venturi_stages_per_path = topology
        .venturi_placements
        .first()
        .map_or(0_i32, |placement| {
            i32::from(placement.serial_throat_count.max(1))
        });
    let strongest_venturi = candidate
        .evaluation
        .venturi
        .placements
        .iter()
        .min_by(|left, right| left.cavitation_number.total_cmp(&right.cavitation_number));
    let cavitation_number = strongest_venturi?.cavitation_number;
    if !cavitation_number.is_finite() {
        return None;
    }

    let cavitation_strength = cavitation_strength_from_sigma(cavitation_number);
    if serial_venturi_stages_per_path <= 0 {
        Some(
            candidate
                .evaluation
                .residence
                .treatment_flow_fraction
                .clamp(0.0, 1.0),
        )
    } else {
        Some(
            candidate
                .evaluation
                .venturi
                .venturi_flow_fraction
                .clamp(0.0, 1.0)
                * (1.0 - (1.0 - cavitation_strength).powi(serial_venturi_stages_per_path)),
        )
    }
}

fn is_hydrosdt_search_candidate(candidate: &BlueprintRankedCandidate) -> bool {
    let Ok(topology) = candidate.candidate.topology_spec() else {
        return false;
    };
    let serial_venturi_stages_per_path = topology
        .venturi_placements
        .first()
        .map_or(0_u8, |placement| placement.serial_throat_count.max(1));
    let strongest_venturi = candidate
        .evaluation
        .venturi
        .placements
        .iter()
        .min_by(|left, right| left.cavitation_number.total_cmp(&right.cavitation_number));
    let Some(strongest_venturi) = strongest_venturi else {
        return false;
    };

    !topology.venturi_placements.is_empty()
        && serial_venturi_stages_per_path > 0
        && strongest_venturi.cavitation_number.is_finite()
        && strongest_venturi.cavitation_number < 1.0
        && candidate.evaluation.separation.cancer_center_fraction > 0.0
}

fn hydrosdt_selection_penalty(
    candidate: &BlueprintRankedCandidate,
    policy: HydroShortlistPolicy,
) -> f64 {
    const NON_HYDROSDT_PENALTY: f64 = 0.25;
    const NO_CAVITATION_GAIN_PENALTY: f64 = 0.10;

    if !is_hydrosdt_search_candidate(candidate) {
        return NON_HYDROSDT_PENALTY;
    }

    let Some(serial_dose_fraction) = hydrosdt_serial_cavitation_dose_fraction(candidate) else {
        return NON_HYDROSDT_PENALTY;
    };
    let shortfall = policy.baseline_serial_cavitation_dose_fraction + policy.required_gain_margin
        - serial_dose_fraction;
    if shortfall <= 0.0 {
        0.0
    } else {
        NO_CAVITATION_GAIN_PENALTY + shortfall
    }
}

fn geometry_concentration_penalty(candidate: &BlueprintRankedCandidate) -> f64 {
    let topology = candidate
        .candidate
        .topology_spec()
        .expect("ranked GA candidate must retain topology metadata");
    let mut lane_family_counts: HashMap<(String, String), usize> = HashMap::new();
    let mut geometry_event_count = 0_usize;
    if let Some(lineage) = candidate.candidate.blueprint.lineage() {
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
                let repeat_depth = (count.saturating_sub(1)) as f64;
                let ancestry_share = count as f64 / geometry_event_count.max(1) as f64;
                let family_bias = if family == "serpentine_compact" {
                    1.2
                } else {
                    1.0
                };
                repeat_depth
                    * lane_serpentine_diminishing_return(route)
                    * family_bias
                    * (1.0 + ancestry_share * geometry_depth_scale)
            })
        })
        .fold(0.0_f64, f64::max)
}

fn operating_point_diversity_penalty(candidate: &BlueprintRankedCandidate) -> f64 {
    let mut operating_point_family_counts: HashMap<String, usize> = HashMap::new();
    if let Some(lineage) = candidate.candidate.blueprint.lineage() {
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

fn ancestry_adjusted_selection_score(
    candidate: &BlueprintRankedCandidate,
    goal: OptimizationGoal,
    hydrosdt_policy: Option<HydroShortlistPolicy>,
) -> f64 {
    let raw_score = evaluation_score_or_zero(&candidate.evaluation);
    let geometry_penalty = geometry_concentration_penalty(candidate);
    let operating_point_penalty = operating_point_diversity_penalty(candidate);
    let hydrosdt_penalty = match (goal, hydrosdt_policy) {
        (OptimizationGoal::InPlaceDeanSerpentineRefinement, Some(policy)) => {
            hydrosdt_selection_penalty(candidate, policy)
        }
        _ => 0.0,
    };
    raw_score
        - geometry_penalty_weight() * geometry_penalty
        - operating_point_penalty_weight() * operating_point_penalty
        - hydrosdt_penalty
}

impl BlueprintGeneticOptimizer {
    #[must_use]
    pub fn new(goal: OptimizationGoal) -> Self {
        Self {
            goal,
            population: 24,
            max_generations: 6,
            top_k: 5,
            seeds: Vec::new(),
            hydrosdt_baseline_serial_cavitation_dose_fraction: None,
        }
    }

    #[must_use]
    pub fn with_population(mut self, population: usize) -> Self {
        self.population = population.max(1);
        self
    }

    #[must_use]
    pub fn with_max_generations(mut self, max_generations: usize) -> Self {
        self.max_generations = max_generations.max(1);
        self
    }

    #[must_use]
    pub fn with_top_k(mut self, top_k: usize) -> Self {
        self.top_k = top_k.max(1);
        self
    }

    #[must_use]
    pub fn with_seeds(mut self, seeds: Vec<BlueprintCandidate>) -> Self {
        self.seeds = seeds;
        self
    }

    #[must_use]
    pub fn with_hydrosdt_baseline_serial_cavitation_dose_fraction(
        mut self,
        baseline_serial_cavitation_dose_fraction: f64,
    ) -> Self {
        self.hydrosdt_baseline_serial_cavitation_dose_fraction =
            Some(baseline_serial_cavitation_dose_fraction.max(0.0));
        self
    }

    pub fn run(&self) -> Result<BlueprintEvolutionResult, OptimError> {
        let mut population = self.initial_population()?;
        let mut archive: HashMap<CandidateFingerprint, BlueprintRankedCandidate> = HashMap::new();
        let mut best_per_generation = Vec::with_capacity(self.max_generations);
        let hydrosdt_policy = match self.goal {
            OptimizationGoal::InPlaceDeanSerpentineRefinement => {
                self.hydrosdt_baseline_serial_cavitation_dose_fraction.map(
                    |baseline_serial_cavitation_dose_fraction| HydroShortlistPolicy {
                        baseline_serial_cavitation_dose_fraction,
                        required_gain_margin: hydrosdt_cavitation_gain_margin(),
                    },
                )
            }
            _ => None,
        };
        let baseline_scores = self
            .seeds
            .iter()
            .map(|candidate| {
                evaluate_goal(candidate, self.goal)
                    .map(|evaluation| evaluation.score.unwrap_or(0.0))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut eval_cache: HashMap<CandidateFingerprint, BlueprintObjectiveEvaluation> =
            HashMap::with_capacity(self.population * 2);

        // ── Phase 1: Broad exploration (topology mutations) ──
        // Standard GA loop: elites + topology mutations (serpentine, venturi,
        // split-merge). Budget: first 2/3 of generations.
        let explore_gens = (self.max_generations * 2 / 3).max(1);
        let refine_gens = self.max_generations - explore_gens;

        for generation in 0..explore_gens {
            if eval_cache.len() > self.population * 4 {
                eval_cache.clear();
            }

            let ranked = rank_population(self.goal, hydrosdt_policy, &population, &mut eval_cache)?;
            if ranked.is_empty() {
                return Err(OptimError::EmptyCandidates);
            }

            best_per_generation.push(evaluation_score_or_zero(&ranked[0].evaluation));
            retain_archive(&mut archive, ranked.iter().cloned());

            if generation % 20 == 19 && archive.len() > 500 {
                let mut sorted: Vec<(CandidateFingerprint, BlueprintRankedCandidate)> =
                    archive.drain().collect();
                sorted.sort_unstable_by(|(_, a), (_, b)| {
                    b.selection_score
                        .total_cmp(&a.selection_score)
                        .then_with(|| {
                            evaluation_score_or_zero(&b.evaluation)
                                .total_cmp(&evaluation_score_or_zero(&a.evaluation))
                        })
                });
                sorted.truncate(500);
                archive.extend(sorted);
            }

            let elite_count = (self.population / 4).max(1).min(ranked.len());
            let elite_ranked = select_diverse_survivors(&ranked, elite_count);
            let mut next_population: Vec<BlueprintCandidate> = elite_ranked
                .iter()
                .map(|ranked| ranked.candidate.clone())
                .collect();
            let mut seen = population_keys(&next_population);

            let mutation_batches: Vec<Vec<BlueprintCandidate>> = (0..elite_ranked.len())
                .map(|i| generate_ga_mutations(&next_population[i]))
                .collect::<Result<_, _>>()?;
            let mut crossover_batches = Vec::new();
            for left_index in 0..elite_ranked.len() {
                for right_index in (left_index + 1)..elite_ranked.len() {
                    crossover_batches.push(generate_ga_crossover_children(
                        &next_population[left_index],
                        &next_population[right_index],
                    )?);
                }
            }
            let mut mutation_offsets = vec![0_usize; mutation_batches.len()];
            let mut crossover_offsets = vec![0_usize; crossover_batches.len()];
            while next_population.len() < self.population {
                let mut advanced = false;
                advanced |= append_round_robin_batch(
                    &mut next_population,
                    &mut seen,
                    &mutation_batches,
                    &mut mutation_offsets,
                    self.population,
                )?;
                advanced |= append_round_robin_batch(
                    &mut next_population,
                    &mut seen,
                    &crossover_batches,
                    &mut crossover_offsets,
                    self.population,
                )?;
                if !advanced {
                    break;
                }
            }

            if next_population.len() < self.population {
                for ranked_candidate in ranked.iter().skip(elite_count) {
                    let key = candidate_key(&ranked_candidate.candidate)?;
                    if seen.insert(key) {
                        next_population.push(ranked_candidate.candidate.clone());
                    }
                    if next_population.len() >= self.population {
                        break;
                    }
                }
            }

            population = next_population;
        }

        // ── Phase 2: Focused refinement (local operating points + continued topology search) ──
        if refine_gens > 0 {
            let ranked = rank_population(self.goal, hydrosdt_policy, &population, &mut eval_cache)?;
            if !ranked.is_empty() {
                let parent_count = (self.population / 6).clamp(1, ranked.len().min(4));
                let parents = select_diverse_survivors(&ranked, parent_count);
                population = build_refine_population(&parents, self.population)?;
                for generation in 0..refine_gens {
                    let ranked =
                        rank_population(self.goal, hydrosdt_policy, &population, &mut eval_cache)?;
                    if ranked.is_empty() {
                        break;
                    }
                    best_per_generation.push(evaluation_score_or_zero(&ranked[0].evaluation));
                    retain_archive(&mut archive, ranked.iter().cloned());

                    if generation + 1 < refine_gens {
                        let parent_count = (self.population / 6).clamp(1, ranked.len().min(4));
                        let parents = select_diverse_survivors(&ranked, parent_count);
                        population = build_refine_population(&parents, self.population)?;
                    }
                }
            }
        }

        let mut all_candidates: Vec<BlueprintRankedCandidate> = archive.into_values().collect();
        all_candidates.sort_unstable_by(|left, right| {
            right
                .selection_score
                .total_cmp(&left.selection_score)
                .then_with(|| {
                    evaluation_score_or_zero(&right.evaluation)
                        .total_cmp(&evaluation_score_or_zero(&left.evaluation))
                })
                .then_with(|| left.candidate.id.cmp(&right.candidate.id))
        });
        // Bound the retained archive to avoid unbounded memory growth.
        // The GA caller only uses the top ~200 for audit + Pareto data.
        let retain_cap = self.top_k.max(200);
        all_candidates.truncate(retain_cap);
        all_candidates.shrink_to_fit();

        for (index, candidate) in all_candidates.iter_mut().enumerate() {
            candidate.rank = index + 1;
            candidate.evaluation = candidate
                .evaluation
                .clone()
                .with_baseline_scores(&baseline_scores);
        }

        all_candidates.sort_unstable_by(|left, right| {
            right
                .evaluation
                .exceeds_all_baselines
                .unwrap_or(false)
                .cmp(&left.evaluation.exceeds_all_baselines.unwrap_or(false))
                .then_with(|| right.selection_score.total_cmp(&left.selection_score))
                .then_with(|| {
                    evaluation_score_or_zero(&right.evaluation)
                        .total_cmp(&evaluation_score_or_zero(&left.evaluation))
                })
                .then_with(|| left.candidate.id.cmp(&right.candidate.id))
        });

        for (index, candidate) in all_candidates.iter_mut().enumerate() {
            candidate.rank = index + 1;
        }

        if all_candidates.is_empty() {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: 0,
            });
        }

        let top_candidate_count = self.top_k.min(all_candidates.len());

        Ok(BlueprintEvolutionResult {
            top_candidates: all_candidates
                .iter()
                .take(top_candidate_count)
                .cloned()
                .collect(),
            all_candidates,
            best_per_generation,
        })
    }

    fn initial_population(&self) -> Result<Vec<BlueprintCandidate>, OptimError> {
        if self.seeds.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        let mut population = Vec::with_capacity(self.population.max(self.seeds.len()));
        let mut seen = HashSet::new();
        for seed in &self.seeds {
            let key = candidate_key(seed)?;
            if seen.insert(key) {
                population.push(seed.clone());
            }
            if population.len() >= self.population {
                break;
            }
        }

        // Generate seed mutations in round-robin order so one prolific seed
        // cannot monopolize the initial population.
        let mutation_batches: Vec<Vec<BlueprintCandidate>> = population
            .iter()
            .map(generate_ga_mutations)
            .collect::<Result<_, _>>()?;
        let mut mutation_offsets = vec![0_usize; mutation_batches.len()];
        while population.len() < self.population {
            let mut advanced = false;
            for (batch_index, batch) in mutation_batches.iter().enumerate() {
                while mutation_offsets[batch_index] < batch.len() {
                    let mutation = batch[mutation_offsets[batch_index]].clone();
                    mutation_offsets[batch_index] += 1;
                    let key = candidate_key(&mutation)?;
                    if seen.insert(key) {
                        population.push(mutation);
                        advanced = true;
                        break;
                    }
                }
                if population.len() >= self.population {
                    break;
                }
            }
            if !advanced {
                break;
            }
        }

        Ok(population)
    }
}

fn rank_population(
    goal: OptimizationGoal,
    hydrosdt_policy: Option<HydroShortlistPolicy>,
    population: &[BlueprintCandidate],
    eval_cache: &mut HashMap<CandidateFingerprint, BlueprintObjectiveEvaluation>,
) -> Result<Vec<BlueprintRankedCandidate>, OptimError> {
    // Skip candidates that fail evaluation (e.g. duplicate channel geometry
    // from deeply nested GA mutations) instead of aborting the entire run.
    let mut ranked: Vec<BlueprintRankedCandidate> = population
        .iter()
        .filter_map(|candidate| {
            let key = candidate_key(candidate).ok();
            let evaluation = if let Some(cached) = key.and_then(|k| eval_cache.get(&k)) {
                cached.clone()
            } else {
                let eval = evaluate_goal(candidate, goal).ok()?;
                if let Some(k) = key {
                    eval_cache.insert(k, eval.clone());
                }
                eval
            };
            Some(BlueprintRankedCandidate {
                rank: 0,
                candidate: candidate.clone(),
                evaluation,
                selection_score: 0.0,
            })
        })
        .collect();
    ranked.retain(|ranked_candidate| ranked_candidate.evaluation.is_eligible());
    for ranked_candidate in &mut ranked {
        ranked_candidate.selection_score =
            ancestry_adjusted_selection_score(ranked_candidate, goal, hydrosdt_policy);
    }
    ranked.sort_by(|left, right| {
        right
            .selection_score
            .partial_cmp(&left.selection_score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                evaluation_score_or_zero(&right.evaluation)
                    .total_cmp(&evaluation_score_or_zero(&left.evaluation))
            })
            .then_with(|| left.candidate.id.cmp(&right.candidate.id))
    });
    for (index, ranked_candidate) in ranked.iter_mut().enumerate() {
        ranked_candidate.rank = index + 1;
    }
    Ok(ranked)
}

fn retain_archive(
    archive: &mut HashMap<CandidateFingerprint, BlueprintRankedCandidate>,
    ranked: impl IntoIterator<Item = BlueprintRankedCandidate>,
) {
    for candidate in ranked {
        if let Ok(key) = candidate_key(&candidate.candidate) {
            match archive.entry(key) {
                Entry::Vacant(entry) => {
                    entry.insert(candidate);
                }
                Entry::Occupied(mut entry) => {
                    if candidate.selection_score > entry.get().selection_score
                        || (candidate.selection_score == entry.get().selection_score
                            && evaluation_score_or_zero(&candidate.evaluation)
                                > evaluation_score_or_zero(&entry.get().evaluation))
                    {
                        entry.insert(candidate);
                    }
                }
            }
        }
    }
}

fn population_keys(population: &[BlueprintCandidate]) -> HashSet<CandidateFingerprint> {
    population
        .iter()
        .filter_map(|candidate| candidate_key(candidate).ok())
        .collect()
}

fn select_diverse_survivors(
    ranked: &[BlueprintRankedCandidate],
    survivor_count: usize,
) -> Vec<BlueprintRankedCandidate> {
    if ranked.is_empty() || survivor_count == 0 {
        return Vec::new();
    }

    let frontier_len = (survivor_count * 4).max(survivor_count).min(ranked.len());
    let frontier = &ranked[..frontier_len];
    let mut selected = Vec::with_capacity(survivor_count.min(frontier_len));
    let mut used = vec![false; frontier_len];

    selected.push(frontier[0].clone());
    used[0] = true;

    while selected.len() < survivor_count && selected.len() < frontier_len {
        let mut best_index = None;
        let mut best_objective = f64::NEG_INFINITY;
        for (index, candidate) in frontier.iter().enumerate() {
            if used[index] {
                continue;
            }
            let novelty = selected
                .iter()
                .map(|existing| diversity_distance(candidate, existing))
                .fold(f64::INFINITY, f64::min);
            let objective = candidate.selection_score + 0.06 * novelty;
            if objective > best_objective {
                best_objective = objective;
                best_index = Some(index);
            }
        }
        let Some(index) = best_index else {
            break;
        };
        used[index] = true;
        selected.push(frontier[index].clone());
    }

    if selected.len() < survivor_count {
        for (index, candidate) in frontier.iter().enumerate() {
            if used[index] {
                continue;
            }
            selected.push(candidate.clone());
            if selected.len() >= survivor_count {
                break;
            }
        }
    }

    selected
}

fn diversity_signature(candidate: &BlueprintRankedCandidate) -> DiversitySignature {
    let topology = candidate
        .candidate
        .topology_spec()
        .expect("ranked GA candidate must retain topology metadata");
    let throat_width_m = topology
        .venturi_placements
        .iter()
        .map(|placement| placement.throat_geometry.throat_width_m)
        .fold(f64::INFINITY, f64::min);
    let max_dean = candidate
        .evaluation
        .venturi
        .placements
        .iter()
        .map(|placement| placement.dean_number)
        .fold(0.0_f64, f64::max);
    DiversitySignature {
        flow_rate_m3_s: candidate.candidate.operating_point.flow_rate_m3_s,
        inlet_gauge_pa: candidate.candidate.operating_point.inlet_gauge_pa,
        throat_width_m: if throat_width_m.is_finite() {
            throat_width_m
        } else {
            0.0
        },
        serial_venturi_stages: topology.serial_venturi_stages() as f64,
        parallel_venturi_count: topology.parallel_venturi_count() as f64,
        split_stage_count: topology.split_stages.len() as f64,
        serpentine_arm_count: topology.serpentine_arm_count() as f64,
        cavitation_selectivity: candidate.evaluation.venturi.cavitation_selectivity_score,
        separation_efficiency: candidate.evaluation.separation.separation_efficiency,
        healthy_cell_protection: healthy_cell_protection_index(
            candidate.evaluation.venturi.wbc_exposure_fraction,
            1.0 - candidate.evaluation.venturi.rbc_exposure_fraction,
        ),
        residence_time_s: candidate.evaluation.residence.treatment_residence_time_s,
        dean_norm: (max_dean / 100.0).clamp(0.0, 1.0),
    }
}

fn normalized_delta(left: f64, right: f64, scale: f64) -> f64 {
    ((left - right).abs() / scale.max(1.0e-12)).clamp(0.0, 1.0)
}

fn diversity_distance(left: &BlueprintRankedCandidate, right: &BlueprintRankedCandidate) -> f64 {
    let left_sig = diversity_signature(left);
    let right_sig = diversity_signature(right);
    normalized_delta(
        left_sig.flow_rate_m3_s,
        right_sig.flow_rate_m3_s,
        left_sig
            .flow_rate_m3_s
            .abs()
            .max(right_sig.flow_rate_m3_s.abs())
            .max(1.0e-9),
    ) + normalized_delta(
        left_sig.inlet_gauge_pa,
        right_sig.inlet_gauge_pa,
        left_sig
            .inlet_gauge_pa
            .abs()
            .max(right_sig.inlet_gauge_pa.abs())
            .max(1.0),
    ) + normalized_delta(
        left_sig.throat_width_m,
        right_sig.throat_width_m,
        left_sig
            .throat_width_m
            .abs()
            .max(right_sig.throat_width_m.abs())
            .max(1.0e-6),
    ) + normalized_delta(
        left_sig.serial_venturi_stages,
        right_sig.serial_venturi_stages,
        4.0,
    ) + normalized_delta(
        left_sig.parallel_venturi_count,
        right_sig.parallel_venturi_count,
        4.0,
    ) + normalized_delta(left_sig.split_stage_count, right_sig.split_stage_count, 4.0)
        + normalized_delta(
            left_sig.serpentine_arm_count,
            right_sig.serpentine_arm_count,
            4.0,
        )
        + normalized_delta(
            left_sig.cavitation_selectivity,
            right_sig.cavitation_selectivity,
            1.0,
        )
        + normalized_delta(
            left_sig.separation_efficiency,
            right_sig.separation_efficiency,
            1.0,
        )
        + normalized_delta(
            left_sig.healthy_cell_protection,
            right_sig.healthy_cell_protection,
            1.0,
        )
        + normalized_delta(
            left_sig.residence_time_s,
            right_sig.residence_time_s,
            left_sig
                .residence_time_s
                .max(right_sig.residence_time_s)
                .max(1.0),
        )
        + normalized_delta(left_sig.dean_norm, right_sig.dean_norm, 1.0)
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

fn lane_serpentine_diminishing_return(route: &cfd_schematics::ChannelRouteSpec) -> f64 {
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
    let smooth_depth =
        positive_ratio((bend_factor - 1.30) / 0.30) + positive_ratio((length_factor - 1.10) / 0.25);
    let long_depth = positive_ratio((length_factor - 1.30) / 0.30)
        + positive_ratio((segments_factor - 1.25) / 0.75);

    let dominant_family_depth = compact_depth
        .max(dense_depth)
        .max(smooth_depth)
        .max(long_depth)
        .max(0.0);
    (dominant_family_depth - 1.0).max(0.0)
}

fn positive_ratio(value: f64) -> f64 {
    value.max(0.0)
}

fn blueprint_with_operating_point_lineage(
    candidate: &BlueprintCandidate,
    family: String,
    operator: &str,
) -> cfd_schematics::domain::model::NetworkBlueprint {
    let mut blueprint = candidate.blueprint.clone();
    let topology = blueprint.topology.clone();
    let lineage = blueprint.lineage.get_or_insert_with(|| {
        topology
            .as_ref()
            .map_or_else(TopologyLineageMetadata::default, |spec| {
                cfd_schematics::BlueprintTopologyFactory::lineage_for_spec(spec)
            })
    });
    lineage.current_stage =
        cfd_schematics::TopologyOptimizationStage::InPlaceDeanSerpentineRefinement;
    lineage.mutations.push(TopologyLineageEvent {
        stage: cfd_schematics::TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        mutation: format!("family={family};lane=global;operator={operator}"),
        source_blueprint: Some(candidate.id.clone()),
    });
    blueprint
}

fn build_refine_population(
    parents: &[BlueprintRankedCandidate],
    population_limit: usize,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let mut population = Vec::with_capacity(population_limit);
    let mut seen = HashSet::with_capacity(population_limit * 2);

    for parent in parents {
        let key = candidate_key(&parent.candidate)?;
        if seen.insert(key) {
            population.push(parent.candidate.clone());
        }
    }

    let local_batches: Vec<Vec<BlueprintCandidate>> = parents
        .iter()
        .map(|parent| focused_operating_point_perturbations(&parent.candidate, population_limit))
        .collect::<Result<_, _>>()?;
    let topology_batches: Vec<Vec<BlueprintCandidate>> = parents
        .iter()
        .map(|parent| generate_ga_mutations(&parent.candidate))
        .collect::<Result<_, _>>()?;
    let mut crossover_batches = Vec::new();
    for left_index in 0..parents.len() {
        for right_index in (left_index + 1)..parents.len() {
            crossover_batches.push(generate_ga_crossover_children(
                &parents[left_index].candidate,
                &parents[right_index].candidate,
            )?);
        }
    }
    let mut local_offsets = vec![0_usize; local_batches.len()];
    let mut topology_offsets = vec![0_usize; topology_batches.len()];
    let mut crossover_offsets = vec![0_usize; crossover_batches.len()];

    while population.len() < population_limit {
        let mut advanced = false;
        advanced |= append_round_robin_batch(
            &mut population,
            &mut seen,
            &local_batches,
            &mut local_offsets,
            population_limit,
        )?;
        advanced |= append_round_robin_batch(
            &mut population,
            &mut seen,
            &topology_batches,
            &mut topology_offsets,
            population_limit,
        )?;
        advanced |= append_round_robin_batch(
            &mut population,
            &mut seen,
            &crossover_batches,
            &mut crossover_offsets,
            population_limit,
        )?;
        if !advanced {
            break;
        }
    }

    Ok(population)
}

fn focused_operating_point_perturbations(
    seed: &BlueprintCandidate,
    limit: usize,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let q_factors = [0.85, 0.95, 1.05, 1.15];
    let p_factors = [0.85, 0.95, 1.05, 1.15];
    let base_q = seed.operating_point.flow_rate_m3_s;
    let base_p = seed.operating_point.inlet_gauge_pa;
    let base_ht = seed.operating_point.feed_hematocrit;
    let mut perturbations = Vec::with_capacity(limit.min(q_factors.len() * p_factors.len()));
    let mut seen = HashSet::with_capacity(limit.min(q_factors.len() * p_factors.len()));

    for &qf in &q_factors {
        for &pf in &p_factors {
            let blueprint = blueprint_with_operating_point_lineage(
                seed,
                format!(
                    "operating_point_refine_q{:.0}_p{:.0}",
                    qf * 100.0,
                    pf * 100.0
                ),
                "operating_point_refine",
            );
            let candidate = BlueprintCandidate::new(
                format!("{}-refine-q{:.0}-p{:.0}", seed.id, qf * 100.0, pf * 100.0),
                blueprint,
                crate::domain::OperatingPoint {
                    flow_rate_m3_s: base_q * qf,
                    inlet_gauge_pa: base_p * pf,
                    feed_hematocrit: base_ht,
                    patient_context: seed.operating_point.patient_context.clone(),
                },
            );
            let key = candidate_key(&candidate)?;
            if seen.insert(key) {
                perturbations.push(candidate);
            }
            if perturbations.len() >= limit {
                return Ok(perturbations);
            }
        }
    }

    Ok(perturbations)
}

fn append_round_robin_batch(
    population: &mut Vec<BlueprintCandidate>,
    seen: &mut HashSet<CandidateFingerprint>,
    batches: &[Vec<BlueprintCandidate>],
    offsets: &mut [usize],
    limit: usize,
) -> Result<bool, OptimError> {
    let mut advanced = false;
    for (batch_index, batch) in batches.iter().enumerate() {
        while offsets[batch_index] < batch.len() {
            let candidate = batch[offsets[batch_index]].clone();
            offsets[batch_index] += 1;
            let key = candidate_key(&candidate)?;
            if seen.insert(key) {
                population.push(candidate);
                advanced = true;
                break;
            }
        }
        if population.len() >= limit {
            break;
        }
    }
    Ok(advanced)
}

fn candidate_key(candidate: &BlueprintCandidate) -> Result<CandidateFingerprint, OptimError> {
    let topology = serde_json::to_vec(candidate.topology_spec()?)
        .map_err(|error| OptimError::InvalidParameter(error.to_string()))?;
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    topology.hash(&mut hasher);
    Ok(CandidateFingerprint {
        topology_hash: hasher.finish(),
        flow_rate_quantized: (candidate.operating_point.flow_rate_m3_s * 1.0e12).round() as i64,
        inlet_gauge_quantized: (candidate.operating_point.inlet_gauge_pa * 1.0e3).round() as i64,
        hematocrit_quantized: (candidate.operating_point.feed_hematocrit * 1.0e6).round() as i64,
    })
}

#[cfg(test)]
mod tests {
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};

    use super::{generate_ga_mutations, BlueprintGeneticOptimizer};
    use crate::application::search::mutations::seed_option2_candidates;

    #[test]
    fn ga_mutations_preserve_branch_width_conservation_and_roles() {
        let seed = canonical_option2_candidate("seed", operating_point(2.0e-6, 30_000.0, 0.18));
        let seeds = seed_option2_candidates(std::slice::from_ref(&seed));
        assert_eq!(seeds.len(), 1);

        let mutations = generate_ga_mutations(&seed).expect("ga mutations");
        assert!(!mutations.is_empty());
        for mutation in mutations {
            let topology = mutation.topology_spec().expect("topology metadata");
            cfd_schematics::BlueprintTopologyFactory::validate_spec(topology)
                .expect("mutated topology must remain valid");
            assert!(
                !mutation.nodes().is_empty(),
                "nodes should remain directly accessible"
            );
            assert!(
                !mutation.channels().is_empty(),
                "channels should remain directly accessible"
            );
        }
    }

    #[test]
    fn blueprint_genetic_optimizer_returns_blueprint_candidates_with_direct_graph_access() {
        let seed = canonical_option2_candidate("seed", operating_point(2.2e-6, 32_000.0, 0.16));
        let result = BlueprintGeneticOptimizer::new(
            crate::OptimizationGoal::InPlaceDeanSerpentineRefinement,
        )
        .with_population(6)
        .with_max_generations(2)
        .with_top_k(2)
        .with_seeds(vec![seed])
        .run()
        .expect("blueprint genetic optimization should succeed");

        assert_eq!(result.top_candidates.len(), 2);
        assert_eq!(result.best_per_generation.len(), 2);
        for ranked in &result.top_candidates {
            assert!(!ranked.candidate.nodes().is_empty());
            assert!(!ranked.candidate.channels().is_empty());
            assert!(
                ranked.candidate.blueprint().topology_spec().is_some(),
                "topology metadata must stay attached to the blueprint"
            );
        }
    }
}
