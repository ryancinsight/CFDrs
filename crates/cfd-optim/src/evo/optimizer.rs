//! [`GeneticOptimizer`] struct and evolutionary search loop.

use std::collections::HashMap;

use rand::{Rng, SeedableRng};

use crate::{
    error::OptimError,
    metrics::compute_metrics,
    orchestration::RankedDesign,
    scoring::{score_candidate_ga, OptimMode, SdtWeights},
};

use super::{
    decode::{candidate_to_genome, decode_genome},
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
    /// Optional warm-start seeds converted from known-good parametric candidates.
    seeds: Vec<crate::design::DesignCandidate>,
    /// Optional deterministic RNG seed for reproducible GA runs.
    rng_seed: Option<u64>,
}

impl GeneticOptimizer {
    /// Construct a new optimizer for the given mode and weights.
    ///
    /// Default population: **300** individuals (≥ 10 pinned per topology for 27 families).
    /// Default generations: **500** (enables convergence then fine-tuning across all topologies).
    ///
    /// The SBX distribution parameter η_c starts at 2.0 (broad exploration) and
    /// ramps linearly to 15.0 by 75% of generations, then holds for exploitation.
    #[must_use]
    pub fn new(mode: OptimMode, weights: SdtWeights) -> Self {
        Self {
            mode,
            weights,
            pop_size: 300,
            max_generations: 500,
            top_k: 5,
            seeds: Vec::new(),
            rng_seed: None,
        }
    }

    /// Warm-start the GA from a set of known-good [`DesignCandidate`]s.
    ///
    /// Up to `pop_size / 2` seed candidates are encoded into genomes and placed
    /// at the front of the initial population.  Remaining slots are filled with
    /// random individuals (with pinned topology diversity).  This lets the GA
    /// refine parametric-sweep winners rather than searching from scratch.
    ///
    /// # Example
    /// ```ignore
    /// let parametric = build_candidate_space();
    /// // score + filter for feasibles …
    /// let seeds: Vec<_> = top_100_feasible_candidates;
    /// let result = GeneticOptimizer::new(mode, weights).with_seeds(seeds).run()?;
    /// ```
    #[must_use]
    pub fn with_seeds(mut self, seeds: Vec<crate::design::DesignCandidate>) -> Self {
        self.seeds = seeds;
        self
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

    /// Set a deterministic RNG seed for reproducible GA runs.
    #[must_use]
    pub fn with_rng_seed(mut self, seed: u64) -> Self {
        self.rng_seed = Some(seed);
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
        let seed = self.rng_seed.unwrap_or_else(|| rand::random::<u64>());
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

        let n_topos = ALL_EVO_TOPOLOGIES.len();

        // ── Initialise population with guaranteed topology diversity ───────
        let mut population: Vec<MillifluidicGenome> = {
            let mut init = Vec::with_capacity(self.pop_size);

            // Phase A: encode warm-start seeds (up to half the population).
            let n_seeds = self.seeds.len().min(self.pop_size / 2);
            for cand in self.seeds.iter().take(n_seeds) {
                init.push(candidate_to_genome(cand));
            }

            // Phase B: pin one random individual per topology to maintain diversity.
            let n_pinned = n_topos.min(self.pop_size - n_seeds);
            for t in 0..n_pinned {
                let mut genes = [0.0f64; N_GENES];
                for g in &mut genes {
                    *g = rng.gen::<f64>();
                }
                genes[0] = (t as f64 + 0.5) / n_topos as f64;
                init.push(MillifluidicGenome { genes });
            }

            // Phase C: fill remaining slots with fully random individuals.
            while init.len() < self.pop_size {
                let mut genes = [0.0f64; N_GENES];
                for g in &mut genes {
                    *g = rng.gen::<f64>();
                }
                init.push(MillifluidicGenome { genes });
            }
            init
        };

        let p_m = 1.0 / N_GENES as f64;
        // Adaptive SBX η: ramp from broad-exploration (2.0) to fine-tuning (15.0)
        // over the first 75% of generations, then hold at the exploitation value.
        const ETA_SBX_MIN: f64 = 2.0;
        const ETA_SBX_MAX: f64 = 15.0;
        const ETA_MUT: f64 = 10.0;
        const TOURNAMENT_K: usize = 3;

        let mut all_feasible: Vec<(f64, crate::design::DesignCandidate)> = Vec::new();
        let mut topo_hall: HashMap<String, (f64, crate::design::DesignCandidate)> = HashMap::new();
        // Bounded archive capacity: retain only the top candidates to avoid
        // unbounded memory growth across generations.
        let archive_cap = self.top_k * 20;

        let mut best_per_gen: Vec<f64> = Vec::with_capacity(self.max_generations);
        let mut mean_per_gen: Vec<f64> = Vec::with_capacity(self.max_generations);
        let mut feasible_per_gen: Vec<usize> = Vec::with_capacity(self.max_generations);
        let mut topology_diversity_per_gen: Vec<usize> = Vec::with_capacity(self.max_generations);

        // Pre-allocated format buffer to avoid per-candidate String allocation.
        let mut id_buf = String::with_capacity(32);

        // ── Main loop ──────────────────────────────────────────────────────
        for gen in 0..self.max_generations {
            let mut fitness: Vec<f64> = Vec::with_capacity(self.pop_size);
            for (i, genome) in population.iter().enumerate() {
                use std::fmt::Write;
                id_buf.clear();
                let _ = write!(id_buf, "g{gen:03}i{i:03}");
                let cand = decode_genome(genome, &id_buf);
                let s = match compute_metrics(&cand) {
                    Ok(m) => {
                        let score =
                            score_candidate_ga(&m, self.mode, &self.weights, cand.inlet_gauge_pa);
                        if score > 0.0 {
                            let tag = cand.topology.short().to_string();
                            match topo_hall.entry(tag) {
                                std::collections::hash_map::Entry::Vacant(v) => {
                                    v.insert((score, cand.clone()));
                                }
                                std::collections::hash_map::Entry::Occupied(mut o) => {
                                    if score > o.get().0 {
                                        o.insert((score, cand.clone()));
                                    }
                                }
                            }
                            all_feasible.push((score, cand));
                        }
                        score
                    }
                    Err(_) => 0.0,
                };
                fitness.push(s);
            }

            // Prune the feasible archive periodically to bound memory usage.
            if all_feasible.len() > archive_cap * 2 {
                all_feasible.sort_by(|a, b| b.0.total_cmp(&a.0));
                all_feasible.truncate(archive_cap);
            }

            // ── Per-generation convergence tracking ────────────────────────
            let mut n_feasible = 0usize;
            let mut best_this_gen = 0.0_f64;
            let mut sum_feasible = 0.0_f64;
            for &s in &fitness {
                if s > 0.0 {
                    n_feasible += 1;
                    sum_feasible += s;
                    if s > best_this_gen {
                        best_this_gen = s;
                    }
                }
            }
            let mean_this_gen = if n_feasible == 0 {
                0.0
            } else {
                sum_feasible / n_feasible as f64
            };
            best_per_gen.push(best_this_gen);
            mean_per_gen.push(mean_this_gen);
            feasible_per_gen.push(n_feasible);

            let n_distinct_topos = {
                let mut seen = [false; 32]; // 28 topologies fit in a fixed-size array
                for g in &population {
                    let idx = (g.genes[0] * (n_topos as f64 - 1e-9))
                        .floor()
                        .clamp(0.0, (n_topos - 1) as f64) as usize;
                    if idx < seen.len() {
                        seen[idx] = true;
                    }
                }
                seen.iter().filter(|&&b| b).count()
            };
            topology_diversity_per_gen.push(n_distinct_topos);

            // ── Elitism: preserve the top 10% ─────────────────────────────
            let mut indexed: Vec<(usize, f64)> = fitness.iter().copied().enumerate().collect();
            indexed.sort_by(|a, b| b.1.total_cmp(&a.1));

            let n_elites = (self.pop_size / 10).max(1);

            // ── Adaptive η_c: linear ramp during exploration phase ─────────
            let eta_sbx = {
                let explore_end = (self.max_generations * 3) / 4;
                if gen < explore_end {
                    ETA_SBX_MIN + (ETA_SBX_MAX - ETA_SBX_MIN) * (gen as f64 / explore_end as f64)
                } else {
                    ETA_SBX_MAX
                }
            };

            // ── Generate offspring ─────────────────────────────────────────
            let mut offspring: Vec<MillifluidicGenome> = Vec::with_capacity(self.pop_size);
            for &(idx, _) in indexed[..n_elites].iter() {
                offspring.push(population[idx]);
            }

            while offspring.len() < self.pop_size {
                let p1_idx = tournament_select(&fitness, TOURNAMENT_K, &mut rng);
                let p2_idx = tournament_select(&fitness, TOURNAMENT_K, &mut rng);
                let (mut c1_genes, mut c2_genes) = sbx_crossover(
                    &population[p1_idx].genes,
                    &population[p2_idx].genes,
                    eta_sbx,
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
                // Tighter dedup: same topology with nearly identical flow AND width
                // is treated as a duplicate.  Tighter flow threshold (0.1e-8 vs 0.5e-8)
                // preserves more diversity in the top-k output.
                if cand.topology == prev.topology
                    && (cand.flow_rate_m3_s - prev.flow_rate_m3_s).abs() < 0.1e-8
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
