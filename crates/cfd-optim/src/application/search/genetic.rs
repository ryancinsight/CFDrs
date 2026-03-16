use std::collections::{hash_map::Entry, HashMap, HashSet};

use crate::application::objectives::{evaluate_goal, BlueprintObjectiveEvaluation};
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use serde::{Deserialize, Serialize};

pub use super::mutations::{generate_ga_mutations, promote_option1_candidate_to_ga_seed};

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

fn evaluation_score_or_zero(evaluation: &BlueprintObjectiveEvaluation) -> f64 {
    evaluation.score.unwrap_or(0.0)
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
        let baseline_scores = self
            .seeds
            .iter()
            .map(|candidate| {
                evaluate_goal(candidate, self.goal)
                    .map(|evaluation| evaluation.score.unwrap_or(0.0))
            })
            .collect::<Result<Vec<_>, _>>()?;

        let mut eval_cache: HashMap<String, BlueprintObjectiveEvaluation> = HashMap::new();

        for _generation in 0..self.max_generations {
            let ranked = rank_population(self.goal, &population, &mut eval_cache)?;
            if ranked.is_empty() {
                return Err(OptimError::EmptyCandidates);
            }

            best_per_generation.push(evaluation_score_or_zero(&ranked[0].evaluation));
            retain_archive(&mut archive, ranked.iter().cloned());

            // Intermediate truncation: prevent unbounded archive growth.
            // Keep top 500 by score every 20 generations (~160 MB → ~80 MB peak).
            if _generation % 20 == 19 && archive.len() > 500 {
                let mut sorted: Vec<(String, BlueprintRankedCandidate)> =
                    archive.drain().collect();
                sorted.sort_unstable_by(|(_, a), (_, b)| {
                    evaluation_score_or_zero(&b.evaluation)
                        .total_cmp(&evaluation_score_or_zero(&a.evaluation))
                });
                sorted.truncate(500);
                archive.extend(sorted);
            }

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
        all_candidates.sort_unstable_by(|left, right| {
            evaluation_score_or_zero(&right.evaluation)
                .total_cmp(&evaluation_score_or_zero(&left.evaluation))
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

        // Generate mutations from the initial seeds to fill the population.
        // Iterate by index to avoid cloning the entire population vec.
        let base_len = population.len();
        for i in 0..base_len {
            if population.len() >= self.population {
                break;
            }
            for mutation in generate_ga_mutations(&population[i])? {
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

fn rank_population(
    goal: OptimizationGoal,
    population: &[BlueprintCandidate],
    eval_cache: &mut HashMap<String, BlueprintObjectiveEvaluation>,
) -> Result<Vec<BlueprintRankedCandidate>, OptimError> {
    let mut ranked = population
        .iter()
        .map(|candidate| {
            // Re-use cached evaluations for elites carried across generations,
            // avoiding redundant evaluate_goal() calls (~15 ms each).
            let key = candidate_key(candidate).ok();
            let evaluation = if let Some(cached) = key.as_deref().and_then(|k| eval_cache.get(k)) {
                Ok(cached.clone())
            } else {
                let eval = evaluate_goal(candidate, goal)?;
                if let Some(k) = key {
                    eval_cache.insert(k, eval.clone());
                }
                Ok(eval)
            };
            evaluation.map(|evaluation| BlueprintRankedCandidate {
                rank: 0,
                candidate: candidate.clone(),
                evaluation,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    ranked.retain(|ranked_candidate| ranked_candidate.evaluation.is_eligible());
    ranked.sort_by(|left, right| {
        evaluation_score_or_zero(&right.evaluation)
            .partial_cmp(&evaluation_score_or_zero(&left.evaluation))
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
                    if evaluation_score_or_zero(&candidate.evaluation)
                        > evaluation_score_or_zero(&entry.get().evaluation)
                    {
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

    use crate::application::search::mutations::seed_option2_candidates;
    use super::{generate_ga_mutations, BlueprintGeneticOptimizer};

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
