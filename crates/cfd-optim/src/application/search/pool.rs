//! Pre-computed evaluation pool for memory-efficient multi-goal ranking.
//!
//! [`EvaluatedPool`] evaluates all blueprint candidates **once** (the expensive
//! 1-D network solve), then allows scoring under different [`OptimizationGoal`]s
//! without re-evaluation — only the cheap scoring function runs.
//!
//! # Memory optimisation
//!
//! - **Scoring cache**: A contiguous `Vec<ScoringSnapshot>` (96 bytes/entry)
//!   holds all scalar fields needed by every scoring function.  The scoring
//!   pass scans this cache sequentially — no pointer-chasing through heap-
//!   allocated metric sub-structs.
//! - **Argsort pattern**: scoring builds a compact `Vec<(u32, f64)>` (8 bytes
//!   per entry), sorts that, and only materialises full `BlueprintObjectiveEvaluation`
//!   objects for the final top-K results.
//! - **Arc sharing**: full evaluations are wrapped in `Arc<BlueprintEvaluation>`
//!   so that `top_k` and `rank` produce output without deep-cloning the
//!   heap-allocated `Vec` fields inside each evaluation.
//!
//! # Example
//!
//! ```rust,ignore
//! use cfd_optim::{EvaluatedPool, OptimizationGoal};
//!
//! let pool = EvaluatedPool::from_candidates(&candidates);
//! let opt1 = pool.top_k(10, OptimizationGoal::AsymmetricSplitResidenceSeparation);
//! let opt2 = pool.top_k(10, OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity);
//! ```

use std::sync::Arc;

use rayon::prelude::*;

use crate::application::orchestration::ScanProgress;
use crate::application::objectives::BlueprintObjectiveEvaluation;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use crate::metrics::{evaluate_blueprint_candidate, BlueprintEvaluation};

// ---------------------------------------------------------------------------
// ScoringSnapshot — compact, cache-friendly scoring data
// ---------------------------------------------------------------------------

/// All scalar fields needed by any [`OptimizationGoal`] scoring function,
/// packed into a single cache-line-friendly struct (~96 bytes).
///
/// Built once during pool construction and scanned sequentially during the
/// scoring pass.  This eliminates pointer-chasing through four separate
/// heap-allocated metric sub-structs for every candidate.
#[derive(Clone, Copy)]
#[repr(C)]
struct ScoringSnapshot {
    // -- Option 1 (selective acoustic residence separation) ---------------
    treatment_residence_time_s: f64,
    treatment_flow_fraction: f64,
    separation_efficiency: f64,
    cancer_center_fraction: f64,
    wbc_center_fraction: f64,
    rbc_peripheral_fraction: f64,
    main_channel_margin: f64,

    // -- Option 2 (selective venturi cavitation) --------------------------
    cavitation_selectivity_score: f64,
    rbc_exposure_fraction: f64,
    wbc_exposure_fraction: f64,
    cavitation_safety_margin: f64,

    // -- GA refinement extras --------------------------------------------
    max_dean_number: f64,
    lineage_mutation_count: u16,
    has_venturi: bool,
    // 5 bytes implicit padding to 96 bytes total
}

// Compile-time size check: ScoringSnapshot must stay ≤ 2 cache lines (128 bytes).
const _: () = assert!(std::mem::size_of::<ScoringSnapshot>() <= 128);

impl ScoringSnapshot {
    /// Build from one candidate/evaluation pair.
    fn from_candidate_eval(candidate: &BlueprintCandidate, eval: &BlueprintEvaluation) -> Self {
        let max_dean_number = eval
            .venturi
            .placements
            .iter()
            .map(|p| p.dean_number)
            .fold(0.0_f64, f64::max);
        Self {
            treatment_residence_time_s: eval.residence.treatment_residence_time_s,
            treatment_flow_fraction: eval.residence.treatment_flow_fraction,
            separation_efficiency: eval.separation.separation_efficiency,
            cancer_center_fraction: eval.separation.cancer_center_fraction,
            wbc_center_fraction: eval.separation.wbc_center_fraction,
            rbc_peripheral_fraction: eval.separation.rbc_peripheral_fraction,
            main_channel_margin: eval.safety.main_channel_margin,
            cavitation_selectivity_score: eval.venturi.cavitation_selectivity_score,
            rbc_exposure_fraction: eval.venturi.rbc_exposure_fraction,
            wbc_exposure_fraction: eval.venturi.wbc_exposure_fraction,
            cavitation_safety_margin: eval.safety.cavitation_safety_margin,
            max_dean_number,
            lineage_mutation_count: candidate.blueprint.lineage().map_or(0, |lineage| {
                lineage.mutations.len().min(u16::MAX as usize) as u16
            }),
            has_venturi: !eval.venturi.placements.is_empty(),
        }
    }

    /// Score this snapshot under the given goal.
    ///
    /// Returns `None` only for *structural* incompatibility (e.g. a non-venturi
    /// candidate under a venturi-required goal).  Physics-based metrics that
    /// are zero or negative reduce the score but never eliminate the candidate
    /// — this preserves gradient signal for optimizers and matches the behavior
    /// of the objective evaluators in `application/objectives/`.
    #[inline]
    fn score(self, goal: OptimizationGoal) -> Option<f64> {
        match goal {
            OptimizationGoal::AsymmetricSplitResidenceSeparation => self.score_option1(),
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity => self.score_option2(),
            OptimizationGoal::InPlaceDeanSerpentineRefinement => self.score_ga(),
        }
    }

    /// Additive scoring for the asymmetric-split residence-separation goal.
    ///
    /// Mirrors [`evaluate_selective_acoustic_residence_separation`] using only
    /// the compact snapshot fields.  Always produces a score in [0.001, 1.0]
    /// — zero-valued metrics reduce the score through the additive weights
    /// but never eliminate the candidate entirely.
    #[inline]
    fn score_option1(self) -> Option<f64> {
        if self.treatment_flow_fraction <= 0.0
            || self.treatment_residence_time_s <= 0.0
            || self.main_channel_margin <= 0.0
        {
            return None;
        }
        let residence_norm = (self.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
        let flow_frac = self.treatment_flow_fraction.clamp(0.0, 1.0);
        let sep = self.separation_efficiency.clamp(0.0, 1.0);
        let cancer = self.cancer_center_fraction.clamp(0.0, 1.0);
        let wbc_exclusion = (1.0 - self.wbc_center_fraction).clamp(0.0, 1.0);
        let rbc_exclusion = self.rbc_peripheral_fraction.clamp(0.0, 1.0);
        let safety = self.main_channel_margin.clamp(0.0, 1.0);

        let base = 0.22 * cancer
            + 0.18 * sep
            + 0.16 * residence_norm
            + 0.12 * flow_frac
            + 0.12 * wbc_exclusion
            + 0.10 * rbc_exclusion
            + 0.10 * safety;
        let healthy_cell_shielding = (wbc_exclusion * rbc_exclusion).sqrt();
        let synergy = 0.12
            * (sep * cancer * residence_norm.max(0.01) * healthy_cell_shielding.max(0.01))
                .powf(0.25);

        Some((base + synergy).clamp(0.001, 1.0))
    }

    /// Additive scoring for the selective venturi-cavitation goal.
    ///
    /// Mirrors [`evaluate_selective_venturi_cavitation`].  Only filters
    /// structurally invalid non-venturi entries; physics-based metrics that
    /// are zero reduce the score through the weighted sum without eliminating
    /// the candidate — matching the objective evaluator's behavior.
    #[inline]
    fn score_option2(self) -> Option<f64> {
        if !self.has_venturi {
            return None;
        }
        if self.cavitation_safety_margin <= 0.0
            || self.treatment_flow_fraction <= 0.0
            || self.cavitation_selectivity_score <= 0.0
        {
            return None;
        }

        let cav = self.cavitation_selectivity_score.clamp(0.0, 1.0);
        let rbc_shield = (1.0 - self.rbc_exposure_fraction).clamp(0.0, 1.0);
        let wbc_shield = (1.0 - self.wbc_exposure_fraction).clamp(0.0, 1.0);
        let safety = self.cavitation_safety_margin.clamp(0.0, 1.0);
        let sep = self.separation_efficiency.clamp(0.0, 1.0);
        let routing_support = self.treatment_flow_fraction.clamp(0.0, 1.0);

        let base = 0.32 * cav
            + 0.18 * rbc_shield
            + 0.14 * wbc_shield
            + 0.12 * safety
            + 0.08 * sep
            + 0.06 * routing_support;
        let synergy = 0.10 * (cav * rbc_shield * routing_support.max(0.01)).cbrt();

        Some((base + synergy).clamp(0.001, 1.0))
    }

    /// Additive scoring for the in-place Dean serpentine refinement goal (GA).
    ///
    /// Mirrors [`evaluate_blueprint_genetic_refinement`].  Only filters
    /// structurally invalid non-venturi entries; the `main_channel_margin`
    /// contributes through its weight rather than gating to zero — matching
    /// the objective evaluator which always returns a score.
    #[inline]
    fn score_ga(self) -> Option<f64> {
        if !self.has_venturi {
            return None;
        }
        if self.main_channel_margin <= 0.0 {
            return None;
        }

        let cav = self.cavitation_selectivity_score.clamp(0.0, 1.0);
        let flow_frac = self.treatment_flow_fraction.clamp(0.0, 1.0);
        let cancer = self.cancer_center_fraction.clamp(0.0, 1.0);
        let sep = self.separation_efficiency.clamp(0.0, 1.0);
        let residence_norm = (self.treatment_residence_time_s / 1.0).clamp(0.0, 1.0);
        let rbc_shield = (1.0 - self.rbc_exposure_fraction).clamp(0.0, 1.0);
        let wbc_shield = (1.0 - self.wbc_exposure_fraction).clamp(0.0, 1.0);
        let safety = self.main_channel_margin.clamp(0.0, 1.0);
        let dean_norm = (self.max_dean_number / 100.0).clamp(0.0, 1.0);
        let lineage_norm = (f64::from(self.lineage_mutation_count) / 5.0).clamp(0.0, 1.0);
        let acoustic_support = {
            let base = 0.22 * cancer
                + 0.18 * sep
                + 0.16 * residence_norm
                + 0.12 * flow_frac
                + 0.12 * (1.0 - self.wbc_center_fraction).clamp(0.0, 1.0)
                + 0.10 * self.rbc_peripheral_fraction.clamp(0.0, 1.0)
                + 0.10 * safety;
            let healthy_cell_shielding = ((1.0 - self.wbc_center_fraction).clamp(0.0, 1.0)
                * self.rbc_peripheral_fraction.clamp(0.0, 1.0))
                .sqrt();
            let synergy = 0.12
                * (sep
                    * cancer
                    * residence_norm.max(0.01)
                    * healthy_cell_shielding.max(0.01))
                .powf(0.25);
            (base + synergy).clamp(0.0, 1.0)
        };

        let base = 0.34 * acoustic_support
            + 0.16 * cav
            + 0.10 * cancer
            + 0.08 * sep
            + 0.08 * rbc_shield
            + 0.06 * wbc_shield
            + 0.06 * safety
            + 0.07 * dean_norm
            + 0.05 * lineage_norm;
        let synergy = 0.18
            * (acoustic_support.max(0.01)
                * cav.max(0.01)
                * cancer.max(0.01)
                * flow_frac.max(0.01)
                * residence_norm.max(0.01)
                * rbc_shield.max(0.01)
                * dean_norm.max(0.01))
            .powf(0.2);

        Some((base + synergy).clamp(0.001, 1.0))
    }
}

// ---------------------------------------------------------------------------
// EvaluatedPool
// ---------------------------------------------------------------------------

/// Lightweight identity for a candidate — just the two strings that
/// `BlueprintObjectiveEvaluation` needs.  Avoids cloning the full
/// `NetworkBlueprint` (nodes, channels, topology specs, metadata…)
/// which was the primary OOM driver at 80 K candidates.
#[derive(Clone)]
struct CandidateIdentity {
    id: String,
    blueprint_name: String,
}

/// Pre-computed evaluation pool for the full candidate space.
///
/// Stores evaluated entries in two parallel vectors:
/// - **`identities`**: `CandidateIdentity` — only the candidate ID and
///   blueprint name (two `String`s, ~100 bytes each).  The full
///   `BlueprintCandidate` is **not** cloned into the pool; callers keep
///   ownership of the original slice.
/// - **`evaluations`**: `Arc<BlueprintEvaluation>` — physics results.
/// - **`scoring_cache`**: `Vec<ScoringSnapshot>` — compact (~96 bytes),
///   contiguous, `Copy` structs holding every scalar the scoring functions
///   need.  The scoring pass scans this linearly, keeping the CPU
///   prefetcher happy and avoiding pointer-chasing.
pub struct EvaluatedPool {
    identities: Vec<CandidateIdentity>,
    evaluations: Vec<Arc<BlueprintEvaluation>>,
    scoring_cache: Vec<ScoringSnapshot>,
}

impl EvaluatedPool {
    /// Build a pool by evaluating all candidates in parallel.
    ///
    /// Only the candidate ID and blueprint name are copied into the pool —
    /// the full `NetworkBlueprint` (nodes, channels, topology specs, …) is
    /// **not** cloned, saving ~10-20 KB per candidate (800 MB–1.6 GB at
    /// 80 K candidates).
    ///
    /// Candidates that fail physics evaluation are silently dropped.
    /// The pool is shrunk to fit after construction to release excess capacity.
    pub fn from_candidates(candidates: &[BlueprintCandidate]) -> Self {
        Self::build_pool(candidates, None)
    }

    /// Build a pool by evaluating all candidates in parallel with heartbeat logging.
    #[must_use]
    pub fn from_candidates_with_progress(
        candidates: &[BlueprintCandidate],
        progress_label: &'static str,
    ) -> Self {
        Self::build_pool(candidates, Some(progress_label))
    }

    fn build_pool(
        candidates: &[BlueprintCandidate],
        progress_label: Option<&'static str>,
    ) -> Self {
        let progress =
            progress_label.map(|label| Arc::new(ScanProgress::new(label, candidates.len())));
        let raw: Vec<(CandidateIdentity, BlueprintEvaluation, ScoringSnapshot)> = candidates
            .par_iter()
            .filter_map(|c| {
                let eval = evaluate_blueprint_candidate(c).ok()?;
                if let Some(progress) = &progress {
                    progress.record();
                }
                let identity = CandidateIdentity {
                    id: c.id.clone(),
                    blueprint_name: c.blueprint.name.clone(),
                };
                let snapshot = ScoringSnapshot::from_candidate_eval(c, &eval);
                Some((identity, eval, snapshot))
            })
            .collect();

        if let Some(progress) = &progress {
            progress.finish();
        }

        let mut identities = Vec::with_capacity(raw.len());
        let mut evaluations = Vec::with_capacity(raw.len());
        let mut scoring_cache = Vec::with_capacity(raw.len());

        for (identity, eval, snapshot) in raw {
            identities.push(identity);
            evaluations.push(Arc::new(eval));
            scoring_cache.push(snapshot);
        }

        identities.shrink_to_fit();
        evaluations.shrink_to_fit();
        scoring_cache.shrink_to_fit();

        Self {
            identities,
            evaluations,
            scoring_cache,
        }
    }

    /// Number of successfully evaluated entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.identities.len()
    }

    /// Whether the pool is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.identities.is_empty()
    }

    /// Approximate heap bytes used by this pool (identities + evaluations + scoring cache).
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        let id_size = std::mem::size_of::<CandidateIdentity>();
        let eval_size = std::mem::size_of::<Arc<BlueprintEvaluation>>();
        let snapshot_size = std::mem::size_of::<ScoringSnapshot>();
        self.identities.len() * (id_size + eval_size + snapshot_size)
    }

    /// Count entries with positive score under the given goal, without
    /// materialising any output objects.
    #[must_use]
    pub fn count_eligible(&self, goal: OptimizationGoal) -> usize {
        self.scoring_cache
            .iter()
            .filter(|snapshot| snapshot.score(goal).is_some())
            .count()
    }

    /// Score and rank ALL entries under the given goal, returning full
    /// evaluations sorted by score descending.
    ///
    /// Uses the argsort pattern internally: builds a compact index+score
    /// vector, sorts it, then materialises the full output.
    pub fn rank(
        &self,
        goal: OptimizationGoal,
    ) -> Result<Vec<BlueprintObjectiveEvaluation>, OptimError> {
        self.top_k(self.identities.len(), goal)
    }

    /// Score all entries and return the top `k` results for the given goal.
    ///
    /// **Phase 1** scans the contiguous `scoring_cache` (no heap chasing).
    /// **Phase 2** sorts the compact `(u32, f64)` argsort vector.
    /// **Phase 3** materialises only the final `k` results, sharing the
    /// `Arc<BlueprintEvaluation>` rather than deep-cloning.
    pub fn top_k(
        &self,
        k: usize,
        goal: OptimizationGoal,
    ) -> Result<Vec<BlueprintObjectiveEvaluation>, OptimError> {
        if self.scoring_cache.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        // Phase 1: score from contiguous cache — no pointer chasing.
        let mut scored: Vec<(u32, f64)> = self
            .scoring_cache
            .iter()
            .enumerate()
            .filter_map(|(idx, snapshot)| snapshot.score(goal).map(|score| (idx as u32, score)))
            .collect();

        if scored.is_empty() {
            return Err(OptimError::EmptyCandidates);
        }

        // Phase 2: sort descending by score.
        scored.sort_unstable_by(|a, b| b.1.total_cmp(&a.1));

        // Phase 3: materialise only the top k via Arc sharing.
        let take = k.min(scored.len());
        let results: Vec<BlueprintObjectiveEvaluation> = scored[..take]
            .iter()
            .map(|&(idx, score)| {
                let i = idx as usize;
                let identity = &self.identities[i];
                let eval = Arc::clone(&self.evaluations[i]);
                BlueprintObjectiveEvaluation::from_identity(
                    goal,
                    &identity.id,
                    &identity.blueprint_name,
                    eval,
                    score,
                )
            })
            .collect();

        Ok(results)
    }
}
