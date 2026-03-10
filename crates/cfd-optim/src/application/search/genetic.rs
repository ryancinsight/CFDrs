use std::collections::{hash_map::Entry, HashMap, HashSet};

use crate::application::objectives::{evaluate_goal, BlueprintObjectiveEvaluation};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use cfd_schematics::{
    topology::presets::with_venturi, BlueprintTopologyFactory, BlueprintTopologyMutation,
    SerpentineSpec, ThroatGeometrySpec, TopologyOptimizationStage,
    VenturiPlacementMode,
};
use cfd_schematics::topology::{TreatmentActuationMode, VenturiConfig};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BlueprintRankedCandidate {
    pub rank: usize,
    pub candidate: BlueprintCandidate,
    pub evaluation: BlueprintObjectiveEvaluation,
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

    pub fn run(&self) -> Result<BlueprintEvolutionResult, OptimError> {
        let mut population = self.initial_population()?;
        let mut archive: HashMap<String, BlueprintRankedCandidate> = HashMap::new();
        let mut best_per_generation = Vec::with_capacity(self.max_generations);

        for _generation in 0..self.max_generations {
            let ranked = rank_population(self.goal, &population)?;
            if ranked.is_empty() {
                return Err(OptimError::EmptyCandidates);
            }

            best_per_generation.push(ranked[0].evaluation.score);
            retain_archive(&mut archive, ranked.iter().cloned());

            let elite_count = (self.population / 4).max(1).min(ranked.len());
            let elites: Vec<BlueprintCandidate> = ranked
                .iter()
                .take(elite_count)
                .map(|ranked| ranked.candidate.clone())
                .collect();

            let mut next_population = elites.clone();
            let mut seen = population_keys(&next_population);

            for elite in &elites {
                for mutation in generate_ga_mutations(elite)? {
                    let key = candidate_key(&mutation)?;
                    if seen.insert(key) {
                        next_population.push(mutation);
                    }
                    if next_population.len() >= self.population {
                        break;
                    }
                }
                if next_population.len() >= self.population {
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

        let mut all_candidates: Vec<BlueprintRankedCandidate> = archive.into_values().collect();
        all_candidates.sort_by(|left, right| {
            right
                .evaluation
                .score
                .partial_cmp(&left.evaluation.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| left.candidate.id.cmp(&right.candidate.id))
        });
        for (index, candidate) in all_candidates.iter_mut().enumerate() {
            candidate.rank = index + 1;
        }

        if all_candidates.len() < self.top_k {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: all_candidates.len(),
            });
        }

        Ok(BlueprintEvolutionResult {
            top_candidates: all_candidates.iter().take(self.top_k).cloned().collect(),
            all_candidates,
            best_per_generation,
        })
    }

    fn initial_population(&self) -> Result<Vec<BlueprintCandidate>, OptimError> {
        let seed_pool = self.seeds.clone();

        if seed_pool.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        let mut population = Vec::with_capacity(self.population.max(seed_pool.len()));
        let mut seen = HashSet::new();
        for seed in seed_pool {
            let has_venturi = seed
                .topology_spec()
                .is_ok_and(|topology| !topology.venturi_placements.is_empty());
            if self.goal == OptimizationGoal::InPlaceDeanSerpentineRefinement && !has_venturi {
                for bootstrap in generate_ga_mutations(&seed)? {
                    let key = candidate_key(&bootstrap)?;
                    if seen.insert(key) {
                        population.push(bootstrap);
                    }
                    if population.len() >= self.population {
                        break;
                    }
                }
            } else {
                let key = candidate_key(&seed)?;
                if seen.insert(key) {
                    population.push(seed);
                }
            }
            if population.len() >= self.population {
                break;
            }
        }

        let base_population = population.clone();
        for seed in &base_population {
            if population.len() >= self.population {
                break;
            }
            for mutation in generate_ga_mutations(seed)? {
                let key = candidate_key(&mutation)?;
                if seen.insert(key) {
                    population.push(mutation);
                }
                if population.len() >= self.population {
                    break;
                }
            }
        }

        Ok(population)
    }
}

pub fn seed_option2_candidates(candidates: &[BlueprintCandidate]) -> Vec<&BlueprintCandidate> {
    candidates
        .iter()
        .filter(|candidate| {
            candidate.inferred_goal()
                == OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity
                || candidate
                    .blueprint
                    .topology_spec()
                    .is_some_and(|spec| !spec.venturi_placements.is_empty())
        })
        .collect()
}

pub fn generate_ga_mutations(
    seed: &BlueprintCandidate,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let topology = seed.topology_spec()?;
    if topology.venturi_placements.is_empty() {
        return bootstrap_ga_venturi_mutations(seed);
    }
    let mut mutated = Vec::new();

    for stage in &topology.split_stages {
        for branch in &stage.branches {
            let width_scale = if branch.treatment_path { 1.08 } else { 0.94 };
            let widened = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::UpdateBranchWidth {
                    stage_id: stage.stage_id.clone(),
                    branch_label: branch.label.clone(),
                    new_width_m: branch.route.width_m * width_scale,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-w-{}", seed.id, branch.label),
                widened,
                seed.operating_point.clone(),
            ));

            if branch.treatment_path {
                let serpentine = branch.route.serpentine.clone().map_or(
                    Some(SerpentineSpec {
                        segments: 4,
                        bend_radius_m: branch.route.width_m * 2.5,
                        segment_length_m: branch.route.length_m / 4.0,
                    }),
                    |serpentine| {
                        Some(SerpentineSpec {
                            segments: serpentine.segments + 1,
                            bend_radius_m: serpentine.bend_radius_m,
                            segment_length_m: serpentine.segment_length_m,
                        })
                    },
                );
                let serpentine_mutation = BlueprintTopologyFactory::mutate(
                    &seed.blueprint,
                    BlueprintTopologyMutation::SetTreatmentSerpentine {
                        stage_id: stage.stage_id.clone(),
                        branch_label: branch.label.clone(),
                        serpentine,
                    },
                    TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                )
                .map_err(OptimError::InvalidParameter)?;
                mutated.push(BlueprintCandidate::new(
                    format!("{}-ga-s-{}", seed.id, branch.label),
                    serpentine_mutation,
                    seed.operating_point.clone(),
                ));
            }
        }
    }

    for placement in &topology.venturi_placements {
        let mut placements = topology.venturi_placements.clone();
        if let Some(target) = placements
            .iter_mut()
            .find(|candidate_placement| candidate_placement.placement_id == placement.placement_id)
        {
            target.throat_geometry.throat_width_m *= 0.92;
            let placement_mutation = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::UpdateVenturiConfiguration {
                    placements,
                    treatment_mode: topology.treatment_mode,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-v-{}", seed.id, placement.placement_id),
                placement_mutation,
                seed.operating_point.clone(),
            ));
        }
    }

    Ok(mutated)
}

fn bootstrap_ga_venturi_mutations(
    seed: &BlueprintCandidate,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let topology = seed.topology_spec()?;
    let treatment_channel_ids = topology.treatment_channel_ids();
    if treatment_channel_ids.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement requires at least one treatment channel, but candidate '{}' has none",
            seed.id
        )));
    }

    let representative_route = topology.channel_route(&treatment_channel_ids[0]).ok_or_else(|| {
        OptimError::InvalidParameter(format!(
            "GA refinement could not resolve treatment channel '{}' in candidate '{}'",
            treatment_channel_ids[0], seed.id
        ))
    })?;
    let placement_mode = if treatment_channel_ids.iter().any(|channel_id| {
        topology
            .channel_route(channel_id)
            .is_some_and(|route| route.serpentine.is_some())
    }) {
        VenturiPlacementMode::CurvaturePeakDeanNumber
    } else {
        VenturiPlacementMode::StraightSegment
    };
    let throat_width_m =
        (representative_route.width_m * 0.4).clamp(60.0e-6, representative_route.width_m * 0.85);
    let throat_length_m = (representative_route.length_m / 8.0)
        .clamp(300.0e-6, representative_route.length_m * 0.5);

    let build_mutant =
        |target_channel_ids: Vec<String>, suffix: &str| -> Result<BlueprintCandidate, OptimError> {
            let augmented = with_venturi(
                topology.clone(),
                VenturiConfig {
                    target_channel_ids,
                    serial_throat_count: 1,
                    throat_geometry: ThroatGeometrySpec {
                        throat_width_m,
                        throat_height_m: representative_route.height_m,
                        throat_length_m,
                        inlet_width_m: representative_route.width_m,
                        outlet_width_m: representative_route.width_m,
                        convergent_half_angle_deg: 7.0,
                        divergent_half_angle_deg: 7.0,
                    },
                    placement_mode,
                },
            );
            let blueprint = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::UpdateVenturiConfiguration {
                    placements: augmented.venturi_placements,
                    treatment_mode: TreatmentActuationMode::VenturiCavitation,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            Ok(BlueprintCandidate::new(
                format!("{}-{suffix}", seed.id),
                blueprint,
                seed.operating_point.clone(),
            ))
        };

    let mut mutated = Vec::with_capacity(treatment_channel_ids.len() + 1);
    mutated.push(build_mutant(
        treatment_channel_ids.clone(),
        "ga-bootstrap-v-all",
    )?);
    for channel_id in treatment_channel_ids {
        mutated.push(build_mutant(
            vec![channel_id.clone()],
            &format!("ga-bootstrap-v-{channel_id}"),
        )?);
    }

    Ok(mutated)
}

fn rank_population(
    goal: OptimizationGoal,
    population: &[BlueprintCandidate],
) -> Result<Vec<BlueprintRankedCandidate>, OptimError> {
    let mut ranked = population
        .iter()
        .map(|candidate| {
            evaluate_goal(candidate, goal).map(|evaluation| BlueprintRankedCandidate {
                rank: 0,
                candidate: candidate.clone(),
                evaluation,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    ranked.sort_by(|left, right| {
        right
            .evaluation
            .score
            .partial_cmp(&left.evaluation.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.candidate.id.cmp(&right.candidate.id))
    });
    for (index, ranked_candidate) in ranked.iter_mut().enumerate() {
        ranked_candidate.rank = index + 1;
    }
    Ok(ranked)
}

fn retain_archive(
    archive: &mut HashMap<String, BlueprintRankedCandidate>,
    ranked: impl IntoIterator<Item = BlueprintRankedCandidate>,
) {
    for candidate in ranked {
        if let Ok(key) = candidate_key(&candidate.candidate) {
            match archive.entry(key) {
                Entry::Vacant(entry) => {
                    entry.insert(candidate);
                }
                Entry::Occupied(mut entry) => {
                    if candidate.evaluation.score > entry.get().evaluation.score {
                        entry.insert(candidate);
                    }
                }
            }
        }
    }
}

fn population_keys(population: &[BlueprintCandidate]) -> HashSet<String> {
    population
        .iter()
        .filter_map(|candidate| candidate_key(candidate).ok())
        .collect()
}

fn candidate_key(candidate: &BlueprintCandidate) -> Result<String, OptimError> {
    let topology = serde_json::to_string(candidate.topology_spec()?)
        .map_err(|error| OptimError::InvalidParameter(error.to_string()))?;
    Ok(format!(
        "{}|{:.12}|{:.3}|{:.6}",
        topology,
        candidate.operating_point.flow_rate_m3_s,
        candidate.operating_point.inlet_gauge_pa,
        candidate.operating_point.feed_hematocrit,
    ))
}

#[cfg(test)]
mod tests {
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};

    use super::{generate_ga_mutations, seed_option2_candidates, BlueprintGeneticOptimizer};

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
