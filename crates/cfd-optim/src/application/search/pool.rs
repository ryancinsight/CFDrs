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
//! let opt1 = pool.top_k(10, OptimizationGoal::SelectiveAcousticResidenceSeparation);
//! let opt2 = pool.top_k(10, OptimizationGoal::SelectiveVenturiCavitation);
//! ```

use std::sync::Arc;

use rayon::prelude::*;

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
    fn new(candidate: &BlueprintCandidate, eval: &BlueprintEvaluation) -> Self {
        let max_dean_number = eval
            .venturi
            .placements
            .iter()
            .map(|p| p.dean_number)
            .fold(0.0_f64, f64::max);
        let lineage_mutation_count = candidate
            .blueprint
            .lineage()
            .map_or(0, |l| l.mutations.len().min(u16::MAX as usize) as u16);
        Self {
            treatment_residence_time_s: eval.residence.treatment_residence_time_s,
            treatment_flow_fraction: eval.residence.treatment_flow_fraction,
            separation_efficiency: eval.separation.separation_efficiency,
            cancer_center_fraction: eval.separation.cancer_center_fraction,
            main_channel_margin: eval.safety.main_channel_margin,
            cavitation_selectivity_score: eval.venturi.cavitation_selectivity_score,
            rbc_exposure_fraction: eval.venturi.rbc_exposure_fraction,
            wbc_exposure_fraction: eval.venturi.wbc_exposure_fraction,
            cavitation_safety_margin: eval.safety.cavitation_safety_margin,
            max_dean_number,
            lineage_mutation_count,
            has_venturi: !eval.venturi.placements.is_empty(),
        }
    }

    /// Score this snapshot under the given goal.
    #[inline]
    fn score(self, goal: OptimizationGoal) -> f64 {
        match goal {
            OptimizationGoal::SelectiveAcousticResidenceSeparation => self.score_option1(),
            OptimizationGoal::SelectiveVenturiCavitation => self.score_option2(),
            OptimizationGoal::BlueprintGeneticRefinement => self.score_ga(),
        }
    }

    #[inline]
    fn score_option1(self) -> f64 {
        let residence =
            self.treatment_residence_time_s * self.treatment_flow_fraction.max(1.0e-9);
        let separation =
            self.separation_efficiency * self.cancer_center_fraction.max(1.0e-9);
        let safety = self.main_channel_margin.max(0.0);
        residence * separation * safety
    }

    #[inline]
    fn score_option2(self) -> f64 {
        if !self.has_venturi {
            return 0.0;
        }
        let cavitation = self.cavitation_selectivity_score.max(0.0);
        let exposure =
            1.0 - 0.5 * (self.rbc_exposure_fraction + self.wbc_exposure_fraction);
        let safety = self.cavitation_safety_margin.max(0.0);
        cavitation * exposure.max(0.0) * safety
    }

    #[inline]
    fn score_ga(self) -> f64 {
        if !self.has_venturi {
            return 0.0;
        }
        let lineage_bonus = self.lineage_mutation_count as f64 * 0.05;
        let dean_bonus = self.max_dean_number / 100.0;
        self.cavitation_selectivity_score.max(0.0)
            * self.separation_efficiency.max(0.0)
            * self.main_channel_margin.max(0.0)
            + lineage_bonus
            + dean_bonus
    }
}

// ---------------------------------------------------------------------------
// EvaluatedPool
// ---------------------------------------------------------------------------

/// Pre-computed evaluation pool for the full candidate space.
///
/// Stores evaluated entries in two parallel vectors:
/// - **`entries`**: `(BlueprintCandidate, Arc<BlueprintEvaluation>)` — full
///   data needed for materialising output results.  `Arc` sharing means
///   `top_k` / `rank` produce results via cheap reference-count bumps
///   rather than deep-cloning heap-allocated Vecs and Strings inside each
///   evaluation.
/// - **`scoring_cache`**: `Vec<ScoringSnapshot>` — compact (~96 bytes),
///   contiguous, `Copy` structs holding every scalar the scoring functions
///   need.  The scoring pass scans this linearly, keeping the CPU
///   prefetcher happy and avoiding pointer-chasing.
pub struct EvaluatedPool {
    entries: Vec<(BlueprintCandidate, Arc<BlueprintEvaluation>)>,
    scoring_cache: Vec<ScoringSnapshot>,
}

impl EvaluatedPool {
    /// Build a pool by evaluating all candidates in parallel.
    ///
    /// Candidates that fail physics evaluation are silently dropped.
    /// The pool is shrunk to fit after construction to release excess capacity.
    pub fn from_candidates(candidates: &[BlueprintCandidate]) -> Self {
        let raw: Vec<(BlueprintCandidate, BlueprintEvaluation)> = candidates
            .par_iter()
            .filter_map(|c| {
                let eval = evaluate_blueprint_candidate(c).ok()?;
                Some((c.clone(), eval))
            })
            .collect();

        let mut entries = Vec::with_capacity(raw.len());
        let mut scoring_cache = Vec::with_capacity(raw.len());

        for (candidate, eval) in raw {
            let snapshot = ScoringSnapshot::new(&candidate, &eval);
            entries.push((candidate, Arc::new(eval)));
            scoring_cache.push(snapshot);
        }

        entries.shrink_to_fit();
        scoring_cache.shrink_to_fit();

        Self {
            entries,
            scoring_cache,
        }
    }

    /// Number of successfully evaluated entries.
    #[must_use]
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Whether the pool is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Approximate heap bytes used by this pool (entries + scoring cache).
    #[must_use]
    pub fn heap_bytes(&self) -> usize {
        let entry_size = std::mem::size_of::<(BlueprintCandidate, Arc<BlueprintEvaluation>)>();
        let snapshot_size = std::mem::size_of::<ScoringSnapshot>();
        self.entries.len() * entry_size + self.scoring_cache.len() * snapshot_size
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
        self.top_k(self.entries.len(), goal)
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
            .filter_map(|(idx, snapshot)| {
                let score = snapshot.score(goal);
                if score > 0.0 {
                    Some((idx as u32, score))
                } else {
                    None
                }
            })
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
                let (candidate, eval) = &self.entries[idx as usize];
                BlueprintObjectiveEvaluation::from_shared_evaluation(
                    goal,
                    candidate,
                    Arc::clone(eval),
                    score,
                )
            })
            .collect();

        Ok(results)
    }
}
