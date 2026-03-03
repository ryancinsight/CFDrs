//! [`GeneticOptimizer`] struct and evolutionary search loop.

use std::collections::HashMap;

use rand::Rng;

use crate::{
    error::OptimError,
    metrics::compute_metrics,
    optimizer::RankedDesign,
    scoring::{score_candidate_ga, OptimMode, SdtWeights},
};

use super::{
    decode::decode_genome,
    operators::{apply_topology_jump, polynomial_mutation, sbx_crossover, tournament_select},
    MillifluidicGenome, ALL_EVO_TOPOLOGIES, N_GENES,
};

/// Full output from a [`GeneticOptimizer`] run, including per-generation
/// convergence statistics for post-run analysis.
#[derive(Debug, Clone)]
pub struct EvolutionResult {
    /// Final top-k ranked designs (diverse, sorted by score descending).
    pub top_designs: Vec<RankedDesign>,
    /// Best feasible score found at each generation (0.0 if none yet feasible).
    pub best_per_gen: Vec<f64>,
    /// Mean feasible score across all feasible designs in each generation.
    pub mean_per_gen: Vec<f64>,
    /// Number of feasible candidates evaluated in each generation.
    pub feasible_per_gen: Vec<usize>,
    /// Number of distinct topology indices present in the population at each
    /// generation.
    pub topology_diversity_per_gen: Vec<usize>,
}

/// Evolutionary optimizer for millifluidic SDT designs covering all topology families.
pub struct GeneticOptimizer {
    mode: OptimMode,
    weights: SdtWeights,
    pop_size: usize,
    max_generations: usize,
    top_k: usize,
}

impl GeneticOptimizer {
    /// Construct a new optimizer for the given mode and weights.
    ///
    /// Default population: **120** individuals (≥ 1 pinned per topology + extras).
    /// Default generations: **200** (sufficient for convergence across 25 topologies).
    #[must_use]
    pub fn new(mode: OptimMode, weights: SdtWeights) -> Self {
        Self {
            mode,
            weights,
            pop_size: 120,
            max_generations: 200,
            top_k: 5,
        }
    }

    /// Override the population size (default: 120).
    #[must_use]
    pub fn with_population(mut self, pop_size: usize) -> Self {
        self.pop_size = pop_size.max(10);
        self
    }

    /// Override the number of generations (default: 200).
    #[must_use]
    pub fn with_max_generations(mut self, max_gen: usize) -> Self {
        self.max_generations = max_gen.max(1);
        self
    }

    /// Override the number of top designs returned (default: 5).
    #[must_use]
    pub fn with_top_k(mut self, k: usize) -> Self {
        self.top_k = k.max(1);
        self
    }

    /// Run the evolutionary search and return an [`EvolutionResult`].
    ///
    /// Uses smooth-penalty fitness for GA navigation (sigmoid ramp near constraint
    /// boundaries instead of hard 0.0 kill), combined with diversity de-duplication
    /// when selecting the final top-k.
    ///
    /// # Errors
    /// Returns [`OptimError::InsufficientCandidates`] if fewer than `top_k`
    /// feasible designs are found after the full search.
    pub fn run(&self) -> Result<EvolutionResult, OptimError> {
        let mut rng = rand::thread_rng();

        let n_topos = ALL_EVO_TOPOLOGIES.len(); // 25

        // ── Initialise population with guaranteed topology diversity ───────
        let mut population: Vec<MillifluidicGenome> = {
            let n_pinned = n_topos.min(self.pop_size);
            let mut seeds = Vec::with_capacity(self.pop_size);
            for t in 0..n_pinned {
                let mut genes: Vec<f64> = (0..N_GENES).map(|_| rng.gen::<f64>()).collect();
                genes[0] = (t as f64 + 0.5) / n_topos as f64;
                seeds.push(MillifluidicGenome { genes });
            }
            for _ in n_pinned..self.pop_size {
                seeds.push(MillifluidicGenome {
                    genes: (0..N_GENES).map(|_| rng.gen::<f64>()).collect(),
                });
            }
            seeds
        };

        let p_m = 1.0 / N_GENES as f64;
        const ETA_SBX: f64 = 2.0;
        const ETA_MUT: f64 = 10.0;
        const TOURNAMENT_K: usize = 3;

        let mut all_feasible: Vec<(f64, crate::design::DesignCandidate)> = Vec::new();
        let mut topo_hall: HashMap<String, (f64, crate::design::DesignCandidate)> = HashMap::new();

        let mut best_per_gen: Vec<f64> = Vec::with_capacity(self.max_generations);
        let mut mean_per_gen: Vec<f64> = Vec::with_capacity(self.max_generations);
        let mut feasible_per_gen: Vec<usize> = Vec::with_capacity(self.max_generations);
        let mut topology_diversity_per_gen: Vec<usize> = Vec::with_capacity(self.max_generations);

        // ── Main loop ──────────────────────────────────────────────────────
        for gen in 0..self.max_generations {
            let fitness: Vec<f64> = population
                .iter()
                .enumerate()
                .map(|(i, genome)| {
                    let cand = decode_genome(genome, &format!("g{gen:03}i{i:03}"));
                    match compute_metrics(&cand) {
                        Ok(m) => {
                            let s = score_candidate_ga(
                                &m,
                                self.mode,
                                &self.weights,
                                cand.inlet_gauge_pa,
                            );
                            if s > 0.0 {
                                let tag = cand.topology.short().to_string();
                                let entry = topo_hall.entry(tag).or_insert((0.0, cand.clone()));
                                if s > entry.0 {
                                    *entry = (s, cand.clone());
                                }
                                all_feasible.push((s, cand));
                            }
                            s
                        }
                        Err(_) => 0.0,
                    }
                })
                .collect();

            // ── Per-generation convergence tracking ────────────────────────
            let gen_feasible_scores: Vec<f64> =
                fitness.iter().copied().filter(|&s| s > 0.0).collect();
            let n_feasible = gen_feasible_scores.len();
            let best_this_gen = gen_feasible_scores.iter().copied().fold(0.0_f64, f64::max);
            let mean_this_gen = if n_feasible == 0 {
                0.0
            } else {
                gen_feasible_scores.iter().sum::<f64>() / n_feasible as f64
            };
            best_per_gen.push(best_this_gen);
            mean_per_gen.push(mean_this_gen);
            feasible_per_gen.push(n_feasible);

            let n_distinct_topos = {
                let mut topo_set: Vec<usize> = population
                    .iter()
                    .map(|g| {
                        (g.genes[0] * (n_topos as f64 - 1e-9))
                            .floor()
                            .clamp(0.0, (n_topos - 1) as f64) as usize
                    })
                    .collect();
                topo_set.sort_unstable();
                topo_set.dedup();
                topo_set.len()
            };
            topology_diversity_per_gen.push(n_distinct_topos);

            // ── Elitism: preserve the top 10% ─────────────────────────────
            let mut indexed: Vec<(usize, f64)> = fitness.iter().copied().enumerate().collect();
            indexed.sort_by(|a, b| b.1.total_cmp(&a.1));

            let n_elites = (self.pop_size / 10).max(1);
            let elites: Vec<MillifluidicGenome> = indexed[..n_elites]
                .iter()
                .map(|(idx, _)| population[*idx].clone())
                .collect();

            // ── Generate offspring ─────────────────────────────────────────
            let mut offspring: Vec<MillifluidicGenome> = elites.clone();

            while offspring.len() < self.pop_size {
                let p1_idx = tournament_select(&fitness, TOURNAMENT_K, &mut rng);
                let p2_idx = tournament_select(&fitness, TOURNAMENT_K, &mut rng);
                let (mut c1_genes, mut c2_genes) = sbx_crossover(
                    &population[p1_idx].genes,
                    &population[p2_idx].genes,
                    ETA_SBX,
                    &mut rng,
                );
                polynomial_mutation(&mut c1_genes, ETA_MUT, p_m, &mut rng);
                polynomial_mutation(&mut c2_genes, ETA_MUT, p_m, &mut rng);
                apply_topology_jump(&mut c1_genes, &mut rng);
                apply_topology_jump(&mut c2_genes, &mut rng);
                offspring.push(MillifluidicGenome { genes: c1_genes });
                if offspring.len() < self.pop_size {
                    offspring.push(MillifluidicGenome { genes: c2_genes });
                }
            }

            population = offspring;
        }

        // ── Diverse top-k selection ────────────────────────────────────────
        for (_, (score, cand)) in topo_hall {
            all_feasible.push((score, cand));
        }

        all_feasible.sort_by(|a, b| b.0.total_cmp(&a.0));
        all_feasible.dedup_by(|a, b| a.0 == b.0);

        let mut diverse_top: Vec<(f64, crate::design::DesignCandidate)> =
            Vec::with_capacity(self.top_k);
        'outer: for (s, cand) in &all_feasible {
            for (_, prev) in &diverse_top {
                if cand.topology == prev.topology
                    && (cand.flow_rate_m3_s - prev.flow_rate_m3_s).abs() < 0.5e-8
                    && (cand.channel_width_m - prev.channel_width_m).abs() < 0.5e-3
                {
                    continue 'outer;
                }
            }
            diverse_top.push((*s, cand.clone()));
            if diverse_top.len() == self.top_k {
                break;
            }
        }

        if diverse_top.len() < self.top_k {
            for (s, cand) in &all_feasible {
                if diverse_top.len() == self.top_k {
                    break;
                }
                if diverse_top.iter().all(|(_, p)| p.id != cand.id) {
                    diverse_top.push((*s, cand.clone()));
                }
            }
        }

        if diverse_top.len() < self.top_k {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: diverse_top.len(),
            });
        }

        let top_designs: Vec<RankedDesign> = diverse_top[..self.top_k]
            .iter()
            .enumerate()
            .filter_map(|(rank_idx, (score, cand))| {
                compute_metrics(cand).ok().map(|metrics| RankedDesign {
                    rank: rank_idx + 1,
                    candidate: cand.clone(),
                    metrics,
                    score: *score,
                })
            })
            .collect();

        if top_designs.len() < self.top_k {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: top_designs.len(),
            });
        }

        Ok(EvolutionResult {
            top_designs,
            best_per_gen,
            mean_per_gen,
            feasible_per_gen,
            topology_diversity_per_gen,
        })
    }
}
