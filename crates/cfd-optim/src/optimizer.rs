//! Parametric optimizer for millifluidic SDT device designs.
//!
//! [`SdtOptimizer`] iterates over the full candidate space, evaluates physics
//! metrics for each candidate using `cfd-1d` models, applies the chosen
//! [`OptimMode`] score function, and returns the top-ranked designs.
//!
//! # Typical usage
//!
//! ```rust,no_run
//! use cfd_optim::{SdtOptimizer, OptimMode, SdtWeights};
//!
//! let weights = SdtWeights::default();
//!
//! // Top 5 designs for hydrodynamic-cavitation SDT
//! let top5_cav = SdtOptimizer::new(OptimMode::SdtCavitation, weights)
//!     .top_k(5)
//!     .expect("optimisation failed");
//!
//! for (rank, design) in top5_cav.iter().enumerate() {
//!     println!("#{} {}: score={:.3}", rank + 1, design.candidate.id, design.score);
//! }
//! ```

use serde::{Deserialize, Serialize};

use crate::design::{build_candidate_space, DesignCandidate};
use crate::error::OptimError;
use crate::metrics::{compute_metrics, SdtMetrics};
use crate::scoring::{score_candidate, OptimMode, SdtWeights};

// ── Output type ──────────────────────────────────────────────────────────────

/// A fully-evaluated, scored design candidate.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RankedDesign {
    /// Ranking position (1 = best).
    pub rank: usize,
    /// The evaluated candidate parameters.
    pub candidate: DesignCandidate,
    /// Physics-derived metrics.
    pub metrics: SdtMetrics,
    /// Scalar objective score in \[0, 1\].  Higher = better.
    pub score: f64,
}

impl RankedDesign {
    /// Return a compact one-line summary suitable for terminal output.
    pub fn summary(&self) -> String {
        let c = &self.candidate;
        let m = &self.metrics;
        format!(
            "#{rank:2}  {id:<40}  score={score:.4}  σ={sigma:>8.3}  \
             main_τ={main_tau:>7.1} Pa  fda={fda}  HI={hi:.2e}  \
             cov={cov:.0}%  res={res:.2}s  ΔP={dp:.0} Pa",
            rank = self.rank,
            id = c.id,
            score = self.score,
            sigma = if m.cavitation_number.is_finite() {
                m.cavitation_number
            } else {
                f64::NEG_INFINITY
            },
            main_tau = m.max_main_channel_shear_pa,
            fda = if m.fda_main_compliant { "OK " } else { "!!" },
            hi = m.hemolysis_index_per_pass,
            cov = m.well_coverage_fraction * 100.0,
            res = m.mean_residence_time_s,
            dp = m.total_pressure_drop_pa,
        )
    }
}

// ── Optimizer ────────────────────────────────────────────────────────────────

/// Multi-objective optimizer for millifluidic SDT device designs.
///
/// Evaluates every candidate in the parameter-sweep space (built by
/// [`build_candidate_space`]) using 1-D physics models, then ranks by the
/// chosen objective score.
pub struct SdtOptimizer {
    mode: OptimMode,
    weights: SdtWeights,
}

impl SdtOptimizer {
    /// Construct a new optimizer for the given mode and weights.
    pub fn new(mode: OptimMode, weights: SdtWeights) -> Self {
        Self { mode, weights }
    }

    /// Evaluate all candidates and return the `k` highest-scoring designs.
    ///
    /// Candidates that fail physics evaluation are silently skipped (they
    /// would score 0 anyway due to invalid parameters).  Only designs that
    /// pass both hard constraints (pressure feasibility + FDA main-channel
    /// shear) appear in the ranked output.
    ///
    /// # Errors
    /// - [`OptimError::EmptyCandidates`] — no candidates survive evaluation.
    /// - [`OptimError::InsufficientCandidates`] — fewer than `k` candidates
    ///   have a non-zero score.
    pub fn top_k(&self, k: usize) -> Result<Vec<RankedDesign>, OptimError> {
        let candidates = build_candidate_space();

        let mut evaluated: Vec<RankedDesign> = candidates
            .into_iter()
            .filter_map(|candidate| {
                match compute_metrics(&candidate) {
                    Ok(metrics) => {
                        let score = score_candidate(&metrics, self.mode, &self.weights);
                        Some(RankedDesign {
                            rank: 0, // filled in below
                            candidate,
                            metrics,
                            score,
                        })
                    }
                    // Skip candidates with invalid parameters / physics errors
                    Err(_) => None,
                }
            })
            .collect();

        if evaluated.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        // Sort descending by score, then by cavitation_potential as tiebreaker
        evaluated.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then_with(|| {
                    b.metrics
                        .cavitation_potential
                        .partial_cmp(&a.metrics.cavitation_potential)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });

        // Count positive-score candidates
        let positive = evaluated.iter().filter(|d| d.score > 0.0).count();
        if positive < k {
            return Err(OptimError::InsufficientCandidates {
                requested: k,
                available: positive,
            });
        }

        // Assign ranks and return only the top-k scoring designs
        let top: Vec<RankedDesign> = evaluated
            .into_iter()
            .filter(|d| d.score > 0.0)
            .take(k)
            .enumerate()
            .map(|(i, mut d)| {
                d.rank = i + 1;
                d
            })
            .collect();

        Ok(top)
    }

    /// Convenience wrapper: return the top 5 designs.
    pub fn top_5(&self) -> Result<Vec<RankedDesign>, OptimError> {
        self.top_k(5)
    }

    /// Run a full evaluation and return ALL non-zero-scored designs, ranked.
    pub fn all_ranked(&self) -> Result<Vec<RankedDesign>, OptimError> {
        let candidates = build_candidate_space();
        let n = candidates.len();

        let mut evaluated: Vec<RankedDesign> = candidates
            .into_iter()
            .filter_map(|candidate| {
                compute_metrics(&candidate).ok().map(|metrics| {
                    let score = score_candidate(&metrics, self.mode, &self.weights);
                    RankedDesign {
                        rank: 0,
                        candidate,
                        metrics,
                        score,
                    }
                })
            })
            .collect();

        if evaluated.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        evaluated.sort_by(|a, b| {
            b.score
                .partial_cmp(&a.score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let total = evaluated.len();
        let ranked: Vec<RankedDesign> = evaluated
            .into_iter()
            .filter(|d| d.score > 0.0)
            .enumerate()
            .map(|(i, mut d)| {
                d.rank = i + 1;
                d
            })
            .collect();

        eprintln!(
            "[cfd-optim] evaluated {total} / {n} candidates, \
             {} pass hard constraints",
            ranked.len()
        );

        Ok(ranked)
    }
}

// ── Statistics helper ────────────────────────────────────────────────────────

/// Compute summary statistics over a slice of ranked designs.
#[derive(Debug, Clone)]
pub struct OptimStats {
    /// Total candidates evaluated (including failures).
    pub total_evaluated: usize,
    /// Number with score > 0 (pass hard constraints).
    pub feasible_count: usize,
    /// Fraction of candidates that are feasible.
    pub feasibility_rate: f64,
    /// Score of the best candidate.
    pub best_score: f64,
    /// Score of the worst feasible candidate.
    pub worst_score: f64,
    /// Mean score across feasible candidates.
    pub mean_score: f64,
}

impl OptimStats {
    /// Compute statistics from a slice of ranked designs.
    pub fn from_ranked(ranked: &[RankedDesign], total_evaluated: usize) -> Self {
        let n = ranked.len();
        if n == 0 {
            return Self {
                total_evaluated,
                feasible_count: 0,
                feasibility_rate: 0.0,
                best_score: 0.0,
                worst_score: 0.0,
                mean_score: 0.0,
            };
        }
        let best = ranked.first().map(|d| d.score).unwrap_or(0.0);
        let worst = ranked.last().map(|d| d.score).unwrap_or(0.0);
        let mean = ranked.iter().map(|d| d.score).sum::<f64>() / n as f64;
        Self {
            total_evaluated,
            feasible_count: n,
            feasibility_rate: n as f64 / total_evaluated as f64,
            best_score: best,
            worst_score: worst,
            mean_score: mean,
        }
    }
}
