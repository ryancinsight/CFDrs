//! Parametric optimizer for millifluidic SDT device designs.
//!
//! [`SdtOptimizer`] iterates over the full candidate space, evaluates physics
//! metrics for each candidate using `cfd-1d` models, applies the chosen
//! [`OptimMode`] score function, and returns the top-ranked designs.
//!
//! # Typical usage
//!
//! ```rust,ignore
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
use std::collections::HashMap;

use rayon::prelude::*;

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
    #[must_use]
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

/// Robustness sweep configuration used by
/// [`SdtOptimizer::top_k_robust`](SdtOptimizer::top_k_robust).
///
/// The robust score is computed by evaluating each candidate over the full
/// Cartesian product of these perturbation sets, then blending:
///
/// - mean scenario score (`mean_weight`)
/// - lower-tail quantile score (`1 - mean_weight`)
///
/// A minimum feasible-scenario ratio is enforced to reject brittle designs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobustSweepConfig {
    /// Multipliers applied to `flow_rate_m3_s`.
    pub flow_scales: Vec<f64>,
    /// Multipliers applied to `inlet_gauge_pa`.
    pub gauge_scales: Vec<f64>,
    /// Absolute feed hematocrit values to test.
    pub feed_hematocrit_values: Vec<f64>,
    /// Additive offsets for `trifurcation_center_frac`.
    pub trifurcation_center_offsets: Vec<f64>,
    /// Additive offsets for staged CIF terminal-bifurcation treatment fraction.
    pub cif_bi_treat_offsets: Vec<f64>,
    /// Additive offsets for `trifurcation_left_frac` in asymmetric trifurcation.
    pub trifurcation_left_offsets: Vec<f64>,
    /// Lower-tail quantile used in robust blending (e.g. `0.10` for P10).
    pub lower_quantile: f64,
    /// Weight on the mean score in blended robust score (`0..1`).
    pub mean_weight: f64,
    /// Minimum feasible-scenario ratio required before down-scaling is applied.
    pub min_feasible_ratio: f64,
    /// Bonus for selective shielding in each scenario:
    /// `(1 - rbc_venturi_exposure_fraction) * cancer_dose_fraction`.
    pub selective_bonus_weight: f64,
}

impl Default for RobustSweepConfig {
    fn default() -> Self {
        Self {
            flow_scales: vec![0.85, 1.00, 1.15],
            gauge_scales: vec![0.90, 1.00, 1.10],
            feed_hematocrit_values: vec![0.35, 0.45, 0.55],
            trifurcation_center_offsets: vec![-0.05, 0.00, 0.05],
            cif_bi_treat_offsets: vec![-0.05, 0.00, 0.05],
            trifurcation_left_offsets: vec![-0.05, 0.00, 0.05],
            lower_quantile: 0.10,
            mean_weight: 0.60,
            min_feasible_ratio: 0.70,
            selective_bonus_weight: 0.10,
        }
    }
}

/// Diagnostics for one candidate evaluated with [`RobustSweepConfig`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RobustScoreStats {
    /// Final blended robust score in `[0, 1]`.
    pub robust_score: f64,
    /// Mean score across all stress scenarios.
    pub mean_score: f64,
    /// Lower-tail quantile score (configured by `lower_quantile`).
    pub lower_quantile_score: f64,
    /// Worst-case scenario score.
    pub worst_score: f64,
    /// Best-case scenario score.
    pub best_score: f64,
    /// Fraction of scenarios with positive score (hard-constraint feasible).
    pub feasible_ratio: f64,
    /// Number of stress scenarios evaluated.
    pub scenarios_evaluated: usize,
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

const ROBUST_SHORTLIST_TOP_GLOBAL: usize = 500;
const ROBUST_SHORTLIST_PER_TOPOLOGY: usize = 20;

impl SdtOptimizer {
    /// Construct a new optimizer for the given mode and weights.
    #[must_use]
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
            .into_par_iter()
            .filter_map(|candidate| self.evaluate_positive_design(candidate))
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

        if evaluated.len() < k {
            return Err(OptimError::InsufficientCandidates {
                requested: k,
                available: evaluated.len(),
            });
        }

        // ── Diversity-aware top-k selection ─────────────────────────────────
        // A candidate is accepted into the diverse set only if it differs from
        // all already-accepted candidates in at least one of:
        //   • topology variant
        //   • flow-rate bucket  (±0.5 µL/min = ±0.5e-8 m³/s)
        //   • channel-width bucket (±0.5 mm = ±0.5e-3 m)
        // Falls back to straight score order if fewer than k diverse candidates
        // can be found.
        let mut diverse: Vec<RankedDesign> = Vec::with_capacity(k);
        'outer: for d in evaluated {
            for prev in &diverse {
                if d.candidate.topology == prev.candidate.topology
                    && (d.candidate.flow_rate_m3_s - prev.candidate.flow_rate_m3_s).abs() < 0.5e-8
                    && (d.candidate.channel_width_m - prev.candidate.channel_width_m).abs() < 0.5e-3
                {
                    continue 'outer; // same cluster — skip
                }
            }
            diverse.push(d);
            if diverse.len() == k {
                break;
            }
        }

        if diverse.len() < k {
            return Err(OptimError::InsufficientCandidates {
                requested: k,
                available: diverse.len(),
            });
        }

        // Assign sequential ranks
        for (i, d) in diverse.iter_mut().enumerate() {
            d.rank = i + 1;
        }

        Ok(diverse)
    }

    /// Convenience wrapper: return the top 5 designs.
    pub fn top_5(&self) -> Result<Vec<RankedDesign>, OptimError> {
        self.top_k(5)
    }

    /// Evaluate all candidates under a perturbation sweep and return the `k`
    /// highest robustly-scoring designs.
    ///
    /// This is intended for selective SDT design selection where nominal-only
    /// rankings can be brittle to operating-point shifts (flow, pressure,
    /// hematocrit, and split-ratio tolerance).
    pub fn top_k_robust(
        &self,
        k: usize,
        config: &RobustSweepConfig,
    ) -> Result<Vec<RankedDesign>, OptimError> {
        let candidates = build_candidate_space();
        let mut nominal_ranked: Vec<RankedDesign> = candidates
            .into_par_iter()
            .filter_map(|candidate| self.evaluate_positive_design(candidate))
            .collect();

        if nominal_ranked.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        nominal_ranked.sort_by(|a, b| {
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

        let shortlist = build_robust_shortlist(
            &nominal_ranked,
            ROBUST_SHORTLIST_TOP_GLOBAL,
            ROBUST_SHORTLIST_PER_TOPOLOGY,
        );

        let mut evaluated: Vec<RankedDesign> = shortlist
            .into_par_iter()
            .filter_map(|nominal| {
                let stats = self.robust_score_stats(&nominal.candidate, config)?;
                if stats.robust_score <= 0.0 {
                    return None;
                }
                Some(RankedDesign {
                    rank: 0,
                    candidate: nominal.candidate,
                    metrics: nominal.metrics,
                    score: stats.robust_score,
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
                .then_with(|| {
                    b.metrics
                        .cavitation_potential
                        .partial_cmp(&a.metrics.cavitation_potential)
                        .unwrap_or(std::cmp::Ordering::Equal)
                })
        });

        if evaluated.len() < k {
            return Err(OptimError::InsufficientCandidates {
                requested: k,
                available: evaluated.len(),
            });
        }

        let mut diverse: Vec<RankedDesign> = Vec::with_capacity(k);
        'outer: for d in evaluated {
            for prev in &diverse {
                if d.candidate.topology == prev.candidate.topology
                    && (d.candidate.flow_rate_m3_s - prev.candidate.flow_rate_m3_s).abs() < 0.5e-8
                    && (d.candidate.channel_width_m - prev.candidate.channel_width_m).abs() < 0.5e-3
                {
                    continue 'outer;
                }
            }
            diverse.push(d);
            if diverse.len() == k {
                break;
            }
        }

        if diverse.len() < k {
            return Err(OptimError::InsufficientCandidates {
                requested: k,
                available: diverse.len(),
            });
        }

        for (i, d) in diverse.iter_mut().enumerate() {
            d.rank = i + 1;
        }

        Ok(diverse)
    }

    /// Convenience wrapper: return the top 5 robust designs.
    pub fn top_5_robust(
        &self,
        config: &RobustSweepConfig,
    ) -> Result<Vec<RankedDesign>, OptimError> {
        self.top_k_robust(5, config)
    }

    /// Run a full evaluation and return ALL non-zero-scored designs, ranked.
    pub fn all_ranked(&self) -> Result<Vec<RankedDesign>, OptimError> {
        let candidates = build_candidate_space();
        let n = candidates.len();

        let mut evaluated: Vec<RankedDesign> = candidates
            .into_par_iter()
            .filter_map(|candidate| self.evaluate_positive_design(candidate))
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
        for (i, d) in evaluated.iter_mut().enumerate() {
            d.rank = i + 1;
        }

        eprintln!(
            "[cfd-optim] evaluated {total} / {n} candidates, \
             {} pass hard constraints",
            evaluated.len()
        );

        Ok(evaluated)
    }

    /// Compute robust-score diagnostics for one candidate across a stress sweep.
    ///
    /// Returns `None` only if no scenario points are available (empty config),
    /// otherwise always returns a stats object with score in `[0, 1]`.
    pub fn robust_score_stats(
        &self,
        candidate: &DesignCandidate,
        config: &RobustSweepConfig,
    ) -> Option<RobustScoreStats> {
        let flow_scales_default = [1.0];
        let gauge_scales_default = [1.0];
        let feed_hcts_default = [candidate.feed_hematocrit];
        let tri_offsets_default = [0.0];
        let cif_bi_offsets_default = [0.0];
        let tri_left_offsets_default = [0.0];

        let flow_scales = if config.flow_scales.is_empty() {
            &flow_scales_default[..]
        } else {
            &config.flow_scales
        };
        let gauge_scales = if config.gauge_scales.is_empty() {
            &gauge_scales_default[..]
        } else {
            &config.gauge_scales
        };
        let feed_hcts = if config.feed_hematocrit_values.is_empty() {
            &feed_hcts_default[..]
        } else {
            &config.feed_hematocrit_values
        };
        let tri_offsets = if config.trifurcation_center_offsets.is_empty() {
            &tri_offsets_default[..]
        } else {
            &config.trifurcation_center_offsets
        };
        let cif_bi_offsets = if config.cif_bi_treat_offsets.is_empty() {
            &cif_bi_offsets_default[..]
        } else {
            &config.cif_bi_treat_offsets
        };
        let tri_left_offsets = if config.trifurcation_left_offsets.is_empty() {
            &tri_left_offsets_default[..]
        } else {
            &config.trifurcation_left_offsets
        };

        let mut scenario_scores: Vec<f64> = Vec::new();
        let mut feasible_count = 0usize;

        for &flow_scale in flow_scales {
            for &gauge_scale in gauge_scales {
                for &feed_hct in feed_hcts {
                    for &tri_offset in tri_offsets {
                        for &cif_bi_offset in cif_bi_offsets {
                            for &tri_left_offset in tri_left_offsets {
                                let mut perturbed = candidate.clone();
                                perturbed.flow_rate_m3_s =
                                    (candidate.flow_rate_m3_s * flow_scale.max(1e-9)).max(1e-12);
                                perturbed.inlet_gauge_pa =
                                    (candidate.inlet_gauge_pa * gauge_scale.max(1e-9)).max(1.0);
                                perturbed.feed_hematocrit = feed_hct.clamp(0.01, 0.60);
                                perturbed.trifurcation_center_frac =
                                    (candidate.trifurcation_center_frac + tri_offset)
                                        .clamp(0.25, 0.70);
                                if matches!(
                                    perturbed.topology,
                                    crate::design::DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
                                ) {
                                    perturbed.cif_pretri_center_frac =
                                        (candidate.cif_pretri_center_frac + tri_offset)
                                            .clamp(0.25, 0.70);
                                    perturbed.cif_terminal_tri_center_frac =
                                        (candidate.cif_terminal_tri_center_frac + tri_offset)
                                            .clamp(0.25, 0.70);
                                    perturbed.cif_terminal_bi_treat_frac =
                                        (candidate.cif_terminal_bi_treat_frac + cif_bi_offset)
                                            .clamp(0.50, 0.85);
                                }
                                if matches!(
                                    perturbed.topology,
                                    crate::design::DesignTopology::AsymmetricTrifurcationVenturi
                                ) {
                                    let min_left = 0.08;
                                    let max_left =
                                        (0.85 - perturbed.trifurcation_center_frac).max(min_left);
                                    perturbed.trifurcation_left_frac =
                                        (candidate.trifurcation_left_frac + tri_left_offset)
                                            .clamp(min_left, max_left);
                                }

                                let scenario_score = match compute_metrics(&perturbed) {
                                    Ok(metrics) => {
                                        let base =
                                            score_candidate(&metrics, self.mode, &self.weights);
                                        let selective_bonus = config
                                            .selective_bonus_weight
                                            .clamp(0.0, 1.0)
                                            * metrics.oncology_selectivity_index.clamp(0.0, 1.0);
                                        let total = (base + selective_bonus).clamp(0.0, 1.0);
                                        if total > 0.0 {
                                            feasible_count += 1;
                                        }
                                        total
                                    }
                                    Err(_) => 0.0,
                                };
                                scenario_scores.push(scenario_score);
                            }
                        }
                    }
                }
            }
        }

        if scenario_scores.is_empty() {
            return None;
        }

        scenario_scores.sort_by(f64::total_cmp);
        let scenarios_evaluated = scenario_scores.len();
        let mean_score = scenario_scores.iter().sum::<f64>() / scenarios_evaluated as f64;
        let lower_quantile_score =
            quantile_sorted(&scenario_scores, config.lower_quantile.clamp(0.0, 1.0));
        let worst_score = *scenario_scores.first().unwrap_or(&0.0);
        let best_score = *scenario_scores.last().unwrap_or(&0.0);
        let feasible_ratio = feasible_count as f64 / scenarios_evaluated as f64;

        let mean_weight = config.mean_weight.clamp(0.0, 1.0);
        let raw_blend = mean_weight * mean_score + (1.0 - mean_weight) * lower_quantile_score;
        let min_feasible_ratio = config.min_feasible_ratio.clamp(0.0, 1.0);
        let feasibility_scale = if min_feasible_ratio <= 0.0 {
            1.0
        } else {
            (feasible_ratio / min_feasible_ratio).clamp(0.0, 1.0)
        };
        let robust_score = (raw_blend * feasibility_scale).clamp(0.0, 1.0);

        Some(RobustScoreStats {
            robust_score,
            mean_score,
            lower_quantile_score,
            worst_score,
            best_score,
            feasible_ratio,
            scenarios_evaluated,
        })
    }

    fn evaluate_positive_design(&self, candidate: DesignCandidate) -> Option<RankedDesign> {
        let metrics = compute_metrics(&candidate).ok()?;
        let score = score_candidate(&metrics, self.mode, &self.weights);
        if score <= 0.0 {
            return None;
        }
        Some(RankedDesign {
            rank: 0,
            candidate,
            metrics,
            score,
        })
    }
}

fn build_robust_shortlist(
    nominal_ranked: &[RankedDesign],
    top_global: usize,
    top_per_topology: usize,
) -> Vec<RankedDesign> {
    let mut shortlist: Vec<RankedDesign> = Vec::new();
    let mut seen_ids = std::collections::HashSet::new();

    for item in nominal_ranked.iter().take(top_global) {
        if seen_ids.insert(item.candidate.id.clone()) {
            shortlist.push(item.clone());
        }
    }

    let mut per_topology_counts: HashMap<crate::design::DesignTopology, usize> = HashMap::new();
    for item in nominal_ranked {
        let count = per_topology_counts
            .entry(item.candidate.topology)
            .or_insert(0);
        if *count >= top_per_topology {
            continue;
        }
        if seen_ids.insert(item.candidate.id.clone()) {
            shortlist.push(item.clone());
        }
        *count += 1;
    }

    shortlist
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
    #[must_use]
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
        let best = ranked.first().map_or(0.0, |d| d.score);
        let worst = ranked.last().map_or(0.0, |d| d.score);
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

fn quantile_sorted(sorted_values: &[f64], q: f64) -> f64 {
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
