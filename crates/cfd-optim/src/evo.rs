//! Evolutionary / genetic algorithm optimizer for millifluidic SDT devices.
//!
//! # Overview
//!
//! [`GeneticOptimizer`] searches the design space using a real-coded genetic
//! algorithm with:
//!
//! - **Tournament selection** (tournament size `k = 3`)
//! - **Simulated Binary Crossover (SBX)** (Deb & Agrawal 1995, η = 2)
//! - **Polynomial mutation** (Deb & Deb 2012, η_m = 20, p_m = 1/n_genes)
//!
//! Each individual is a [`MillifluidicGenome`] — a fixed-length real-valued
//! vector that encodes one fully-specified [`DesignCandidate`].
//!
//! The optimizer covers **all 14 fixed topology families** plus the GA-only
//! `AdaptiveTree` (variable-depth split tree), giving **15 searchable families**:
//!
//! | Topology                         | Outlets |
//! |----------------------------------|---------|
//! | `SingleVenturi`                  | 1       |
//! | `BifurcationVenturi`             | 2       |
//! | `TrifurcationVenturi`            | 3       |
//! | `VenturiSerpentine`              | 1       |
//! | `SerpentineGrid`                 | 1       |
//! | `CellSeparationVenturi`          | 2       |
//! | `WbcCancerSeparationVenturi`     | 2       |
//! | `DoubleBifurcationVenturi`       | 4       |
//! | `TripleBifurcationVenturi`       | 8       |
//! | `DoubleTrifurcationVenturi`      | 9       |
//! | `BifurcationTrifurcationVenturi` | 6       |
//! | `SerialDoubleVenturi`            | 1       |
//! | `BifurcationSerpentine`          | 2       |
//! | `TrifurcationSerpentine`         | 3       |
//! | `AdaptiveTree`                   | 1–81    |
//!
//! When gene 0 selects the `AdaptiveTree` slot (index 14), genes 8–12 are
//! decoded to set the tree depth (0–4) and per-level split type (Bi / Tri).
//! The maximum depth is capped so that no leaf channel narrows below 150 µm.
//!
//! # Usage
//!
//! ```rust,no_run
//! use cfd_optim::{OptimMode, SdtWeights, evo::GeneticOptimizer};
//!
//! let best = GeneticOptimizer::new(
//!     OptimMode::SdtCavitation,
//!     SdtWeights::default(),
//! )
//! .with_population(80)
//! .with_max_generations(150)
//! .run()
//! .expect("evolutionary search failed");
//!
//! for d in &best {
//!     println!("{}", d.summary());
//! }
//! ```
//!
//! # References
//! - Deb, K. & Agrawal, R. B. (1995). Simulated binary crossover for continuous
//!   search space. *Complex Systems*, 9, 115–148.
//! - Deb, K. & Deb, D. (2012). Analysing mutation schemes for real-parameter
//!   genetic algorithms. *Int. J. Artif. Intell. Soft Comput.*, 4, 1–28.

use rand::Rng;

use crate::{
    constraints::{
        CHANNEL_HEIGHT_M, FLOW_RATES_M3_S, INLET_GAUGES_PA,
        THROAT_DIAMETERS_M, TREATMENT_WIDTH_MM, VENTURI_INLET_DIAM_M,
    },
    design::{DesignCandidate, DesignTopology},
    error::OptimError,
    metrics::compute_metrics,
    optimizer::RankedDesign,
    scoring::{score_candidate, OptimMode, SdtWeights},
};

// ── All 15 topology families available to the genetic search ─────────────────
// Index 14 is the AdaptiveTree placeholder; actual depth/split_types are
// decoded from genes 8–12 inside `decode_genome` when topo_idx == 14.

const ALL_EVO_TOPOLOGIES: [DesignTopology; 15] = [
    DesignTopology::SingleVenturi,
    DesignTopology::BifurcationVenturi,
    DesignTopology::TrifurcationVenturi,
    DesignTopology::VenturiSerpentine,
    DesignTopology::SerpentineGrid,
    DesignTopology::CellSeparationVenturi,
    DesignTopology::WbcCancerSeparationVenturi,
    DesignTopology::DoubleBifurcationVenturi,
    DesignTopology::TripleBifurcationVenturi,
    DesignTopology::DoubleTrifurcationVenturi,
    DesignTopology::BifurcationTrifurcationVenturi,
    DesignTopology::SerialDoubleVenturi,
    DesignTopology::BifurcationSerpentine,
    DesignTopology::TrifurcationSerpentine,
    // Placeholder — depth and split_types injected from genes 8–12:
    DesignTopology::AdaptiveTree { levels: 0, split_types: 0 },
];

// ── Genome definition ─────────────────────────────────────────────────────────

/// Real-coded genome for a millifluidic design with a single inlet and outlet.
///
/// Gene layout (all normalised to `[0, 1]`):
///
/// | Index | Meaning | Range (physical) |
/// |-------|---------|-----------------|
/// | 0  | Topology index (continuous, rounded to int) | 0 → SingleVenturi … 14 → AdaptiveTree (see `ALL_EVO_TOPOLOGIES`) |
/// | 1  | Flow rate | \[`Q_min`, `Q_max`\] m³/s |
/// | 2  | Inlet gauge pressure | \[`P_min`, `P_max`\] Pa |
/// | 3  | Throat diameter (venturi only) | \[`d_min`, `d_max`\] m |
/// | 4  | Channel width | \[`w_min`, `w_max`\] m |
/// | 5  | Serpentine segment count (serpentine topologies) | 2 – 12 |
/// | 6  | Segment length fraction of treatment width | 0.5 – 1.5 |
/// | 7  | Bend radius fraction of segment length | 0.05 – 0.25 |
/// | 8  | AdaptiveTree depth | 0 – 4 (capped so leaf width ≥ 150 µm) |
/// | 9  | AdaptiveTree level-0 split type | 0 → Bifurcation, 1 → Trifurcation |
/// | 10 | AdaptiveTree level-1 split type | same encoding |
/// | 11 | AdaptiveTree level-2 split type | same encoding |
/// | 12 | AdaptiveTree level-3 split type | same encoding |
///
/// Genes 8–12 are only meaningful when gene 0 selects topology index 14
/// (`AdaptiveTree`).  For all other topologies they are decoded but ignored.
#[derive(Debug, Clone)]
pub struct MillifluidicGenome {
    /// Normalised gene values ∈ [0, 1].
    pub genes: Vec<f64>,
}

/// Number of genes per individual (8 common + 5 AdaptiveTree-specific).
const N_GENES: usize = 13;

/// Minimum leaf channel width for AdaptiveTree depth constraint [m].
///
/// Channels below this width are not fabricatable with standard soft-lithography.
const MIN_ADAPTIVE_CH_M: f64 = 150e-6;

/// Decode a normalised genome into a [`DesignCandidate`].
///
/// Continuous gene 0 is rounded to select the topology.
/// Gene 5 is rounded to select the segment count.
/// When gene 0 selects the `AdaptiveTree` slot (index 14), genes 8–12 are
/// used to set the tree depth and per-level split types; the depth is capped
/// so that no leaf channel falls below [`MIN_ADAPTIVE_CH_M`].
pub fn decode_genome(g: &MillifluidicGenome, id_prefix: &str) -> DesignCandidate {
    let genes = &g.genes;

    // Gene 4: channel width decoded early — needed for AdaptiveTree depth cap.
    // (Linear between 1.0 mm and 4.0 mm — millifluidic range.)
    let w_ch = 1.0e-3 + genes[4] * 3.0e-3;

    // Gene 0: topology
    let topo_idx = (genes[0] * (ALL_EVO_TOPOLOGIES.len() as f64 - 1e-9))
        .floor()
        .clamp(0.0, (ALL_EVO_TOPOLOGIES.len() - 1) as f64) as usize;

    // If the last slot (AdaptiveTree) is selected, decode depth + split_types
    // from genes 8–12.  Otherwise use the fixed topology table directly.
    let topology = if topo_idx == ALL_EVO_TOPOLOGIES.len() - 1 {
        // Compute maximum allowed depth: conservatively use bifurcation fan (2)
        // so that w_ch / 2^max_depth ≥ MIN_ADAPTIVE_CH_M.
        let max_depth = {
            let mut d = 0u8;
            let mut w = w_ch;
            while d < 4 {
                w /= 2.0;                         // worst case: every level is Bi
                if w < MIN_ADAPTIVE_CH_M { break; }
                d += 1;
            }
            d
        };

        // Gene 8: tree depth (0 → max_depth).
        let depth = if max_depth == 0 {
            0u8
        } else {
            ((genes[8] * (max_depth as f64 + 1.0)).floor() as u8).min(max_depth)
        };

        // Genes 9–12: one bit per level (0 = Bi, 1 = Tri).
        let mut split_types = 0u8;
        for i in 0..4usize {
            if genes[9 + i] > 0.5 {
                split_types |= 1u8 << i;
            }
        }
        // Mask off bits for unused levels (depth … 3).
        let mask = if depth == 0 { 0u8 } else { (1u8 << depth).wrapping_sub(1) };
        split_types &= mask;

        DesignTopology::AdaptiveTree { levels: depth, split_types }
    } else {
        ALL_EVO_TOPOLOGIES[topo_idx]
    };

    // Gene 1: flow rate (log-linear between array bounds, constants are sorted ascending)
    let q_min = FLOW_RATES_M3_S[0];
    let q_max = FLOW_RATES_M3_S[FLOW_RATES_M3_S.len() - 1];
    let q = q_min * (q_max / q_min).powf(genes[1]);

    // Gene 2: inlet gauge pressure (linear)
    let p_min = INLET_GAUGES_PA[0];
    let p_max = INLET_GAUGES_PA[INLET_GAUGES_PA.len() - 1];
    let gauge = p_min + genes[2] * (p_max - p_min);

    // Gene 3: throat diameter (linear, only used when topology has a venturi)
    let d_min = THROAT_DIAMETERS_M[0];
    let d_max = THROAT_DIAMETERS_M[THROAT_DIAMETERS_M.len() - 1];
    let d_throat = if topology.has_venturi() {
        d_min + genes[3] * (d_max - d_min)
    } else {
        0.0
    };

    // Gene 5: serpentine segment count (2 – 12, only for serpentine topologies)
    let n_segs = if topology.has_serpentine() {
        2 + (genes[5] * 10.0).round() as usize
    } else {
        1
    };

    // Gene 6: segment length as fraction of treatment width (0.5 – 1.5)
    let seg_len = (0.5 + genes[6]) * TREATMENT_WIDTH_MM * 1e-3;

    // Gene 7: bend radius as fraction of segment length (0.05 – 0.25)
    let bend_r = (0.05 + genes[7] * 0.20) * seg_len;

    let throat_len = if d_throat > 0.0 { d_throat * 2.0 } else { 0.0 };

    // Encode topology identity tag (include depth+split for AdaptiveTree)
    let topo_tag = match topology {
        DesignTopology::AdaptiveTree { levels, split_types } => {
            format!("AT-d{}-s{:04b}", levels, split_types & 0x0F)
        }
        _ => format!("t{}", topo_idx),
    };

    let id = format!(
        "{}-EVO-{}-q{:.0}-g{:.0}-d{:.0}-w{:.0}-n{}",
        id_prefix,
        topo_tag,
        q * 6e7,
        gauge * 1e-3,
        d_throat * 1e6,
        w_ch * 1e6,
        n_segs,
    );

    DesignCandidate {
        id,
        topology,
        flow_rate_m3_s: q,
        inlet_gauge_pa: gauge,
        throat_diameter_m: d_throat,
        inlet_diameter_m: VENTURI_INLET_DIAM_M,
        throat_length_m: throat_len,
        channel_width_m: w_ch,
        channel_height_m: CHANNEL_HEIGHT_M,
        serpentine_segments: n_segs,
        segment_length_m: seg_len,
        bend_radius_m: bend_r,
    }
}

// ── Genetic operators ─────────────────────────────────────────────────────────

/// Tournament selection: return the index of the winner among `k` random participants.
fn tournament_select<R: Rng>(fitnesses: &[f64], k: usize, rng: &mut R) -> usize {
    let n = fitnesses.len();
    let mut best = rng.gen_range(0..n);
    for _ in 1..k {
        let challenger = rng.gen_range(0..n);
        if fitnesses[challenger] > fitnesses[best] {
            best = challenger;
        }
    }
    best
}

/// Simulated Binary Crossover (SBX, Deb & Agrawal 1995).
///
/// Produces two offspring from two parents.  `eta` is the distribution index
/// (larger = children closer to parents).
fn sbx_crossover<R: Rng>(
    p1: &[f64],
    p2: &[f64],
    eta: f64,
    rng: &mut R,
) -> (Vec<f64>, Vec<f64>) {
    let mut c1 = p1.to_vec();
    let mut c2 = p2.to_vec();

    for i in 0..p1.len() {
        if rng.gen::<f64>() < 0.5 {
            let u = rng.gen::<f64>();
            let beta = if u <= 0.5 {
                (2.0 * u).powf(1.0 / (eta + 1.0))
            } else {
                (1.0 / (2.0 * (1.0 - u))).powf(1.0 / (eta + 1.0))
            };
            c1[i] = (0.5 * ((1.0 + beta) * p1[i] + (1.0 - beta) * p2[i])).clamp(0.0, 1.0);
            c2[i] = (0.5 * ((1.0 - beta) * p1[i] + (1.0 + beta) * p2[i])).clamp(0.0, 1.0);
        }
    }
    (c1, c2)
}

/// Polynomial mutation (Deb & Deb 2012).
///
/// Mutates each gene with probability `p_m`.  `eta_m` is the mutation
/// distribution index (larger = smaller perturbations).
fn polynomial_mutation<R: Rng>(genes: &mut Vec<f64>, eta_m: f64, p_m: f64, rng: &mut R) {
    for g in genes.iter_mut() {
        if rng.gen::<f64>() < p_m {
            let u = rng.gen::<f64>();
            let delta = if u < 0.5 {
                (2.0 * u).powf(1.0 / (eta_m + 1.0)) - 1.0
            } else {
                1.0 - (2.0 * (1.0 - u)).powf(1.0 / (eta_m + 1.0))
            };
            *g = (*g + delta).clamp(0.0, 1.0);
        }
    }
}

// ── Optimizer ─────────────────────────────────────────────────────────────────

/// Evolutionary optimizer for millifluidic SDT designs covering all 14 topology families.
///
/// Uses a real-coded genetic algorithm to search beyond the discrete grid of
/// [`build_candidate_space`][crate::design::build_candidate_space].
pub struct GeneticOptimizer {
    mode: OptimMode,
    weights: SdtWeights,
    pop_size: usize,
    max_generations: usize,
    top_k: usize,
}

impl GeneticOptimizer {
    /// Construct a new optimizer for the given mode and weights.
    pub fn new(mode: OptimMode, weights: SdtWeights) -> Self {
        Self {
            mode,
            weights,
            pop_size: 60,
            max_generations: 120,
            top_k: 5,
        }
    }

    /// Override the population size (default: 60).
    pub fn with_population(mut self, pop_size: usize) -> Self {
        self.pop_size = pop_size.max(10);
        self
    }

    /// Override the number of generations (default: 120).
    pub fn with_max_generations(mut self, max_gen: usize) -> Self {
        self.max_generations = max_gen.max(1);
        self
    }

    /// Override the number of top designs returned (default: 5).
    pub fn with_top_k(mut self, k: usize) -> Self {
        self.top_k = k.max(1);
        self
    }

    /// Run the evolutionary search and return the top `k` designs.
    ///
    /// # Errors
    /// Returns [`OptimError::InsufficientCandidates`] if fewer than `top_k`
    /// feasible designs are found after the full search.
    pub fn run(&self) -> Result<Vec<RankedDesign>, OptimError> {
        let mut rng = rand::thread_rng();

        // ── Initialise population ──────────────────────────────────────────
        let mut population: Vec<MillifluidicGenome> = (0..self.pop_size)
            .map(|_| MillifluidicGenome {
                genes: (0..N_GENES).map(|_| rng.gen::<f64>()).collect(),
            })
            .collect();

        let p_m = 1.0 / N_GENES as f64;
        const ETA_SBX: f64 = 2.0;
        const ETA_MUT: f64 = 20.0;
        const TOURNAMENT_K: usize = 3;

        let mut best_scored: Vec<(f64, DesignCandidate)> = Vec::new();

        // ── Main loop ──────────────────────────────────────────────────────
        for gen in 0..self.max_generations {
            // Evaluate fitness for each individual
            let fitness: Vec<f64> = population
                .iter()
                .enumerate()
                .map(|(i, genome)| {
                    let cand = decode_genome(genome, &format!("g{:03}i{:03}", gen, i));
                    match compute_metrics(&cand) {
                        Ok(m) => {
                            let s = score_candidate(&m, self.mode, &self.weights);
                            if s > 0.0 {
                                best_scored.push((s, cand));
                            }
                            s
                        }
                        Err(_) => 0.0,
                    }
                })
                .collect();

            // ── Elitism: preserve the top 10% of current population ───────
            let mut indexed: Vec<(usize, f64)> = fitness
                .iter()
                .copied()
                .enumerate()
                .collect();
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
                offspring.push(MillifluidicGenome { genes: c1_genes });
                if offspring.len() < self.pop_size {
                    offspring.push(MillifluidicGenome { genes: c2_genes });
                }
            }

            population = offspring;
        }

        // ── Collect and rank unique top-k ──────────────────────────────────
        // Sort all feasible candidates found across all generations.
        best_scored.sort_by(|a, b| b.0.total_cmp(&a.0));
        best_scored.dedup_by(|a, b| a.0 == b.0); // remove exact duplicates

        if best_scored.len() < self.top_k {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: best_scored.len(),
            });
        }

        let results: Vec<RankedDesign> = best_scored[..self.top_k]
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

        if results.len() < self.top_k {
            return Err(OptimError::InsufficientCandidates {
                requested: self.top_k,
                available: results.len(),
            });
        }

        Ok(results)
    }
}
