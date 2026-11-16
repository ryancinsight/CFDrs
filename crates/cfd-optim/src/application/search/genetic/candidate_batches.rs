//! Per-parent candidate batches in one contiguous buffer.
//!
//! The genetic search generates a batch of mutations, crossovers, or local
//! perturbations per parent, then appends them round-robin so one prolific
//! parent cannot monopolize the next generation. Those batches were a
//! `Vec<BlueprintCandidate>` per parent: a pointer table plus a heap block per
//! parent, rebuilt every generation, and the round-robin walk existed twice —
//! once as [`append_round_robin_batch`] and once open-coded in the seed path.
//!
//! [`CandidateBatches`] holds every batch back to back behind an offset table,
//! so a generation costs one candidate buffer plus one offset table, and both
//! walks are the one helper below.

use std::collections::HashSet;

use super::{CandidateFingerprint, candidate_key};
use crate::domain::BlueprintCandidate;
use crate::error::OptimError;

/// Candidate batches stored back to back, batch `i` in
/// `candidates[offsets[i]..offsets[i + 1]]`.
#[derive(Debug, Default)]
pub(super) struct CandidateBatches {
    candidates: Vec<BlueprintCandidate>,
    offsets: Vec<usize>,
}

impl CandidateBatches {
    /// An empty set with its leading zero offset.
    pub(super) fn new() -> Self {
        Self {
            candidates: Vec::new(),
            offsets: vec![0],
        }
    }

    /// Assemble from per-parent batch results, appending as each arrives.
    ///
    /// Each batch buffer is consumed and released in turn, so the peak shape
    /// is one flat buffer under growth rather than one heap block per parent.
    pub(super) fn from_results<I>(batches: I) -> Result<Self, OptimError>
    where
        I: IntoIterator<Item = Result<Vec<BlueprintCandidate>, OptimError>>,
    {
        let mut assembled = Self::new();
        for batch in batches {
            assembled.push_batch(batch?);
        }
        Ok(assembled)
    }

    /// Append one batch, recording its end offset.
    pub(super) fn push_batch(&mut self, batch: Vec<BlueprintCandidate>) {
        self.candidates.extend(batch);
        self.offsets.push(self.candidates.len());
    }

    /// The number of batches.
    pub(super) fn batch_count(&self) -> usize {
        self.offsets.len() - 1
    }

    /// Batch `index` in insertion order.
    ///
    /// # Panics
    ///
    /// Panics when `index >= self.batch_count()`.
    pub(super) fn batch(&self, index: usize) -> &[BlueprintCandidate] {
        &self.candidates[self.offsets[index]..self.offsets[index + 1]]
    }

    /// One read cursor per batch, all at the start.
    pub(super) fn cursors(&self) -> Vec<usize> {
        vec![0; self.batch_count()]
    }
}

/// Append one unseen candidate per batch, in batch order, up to `limit`.
///
/// The single round-robin walk: a parent whose first few children duplicate
/// already-seen candidates is skipped forward rather than blocking the rest.
/// Returns whether anything was appended, so the caller can stop once every
/// batch is exhausted.
///
/// # Errors
///
/// Returns [`OptimError`] when a candidate's fingerprint cannot be built.
pub(super) fn append_round_robin_batch(
    population: &mut Vec<BlueprintCandidate>,
    seen: &mut HashSet<CandidateFingerprint>,
    batches: &CandidateBatches,
    cursors: &mut [usize],
    limit: usize,
) -> Result<bool, OptimError> {
    let mut advanced = false;
    for batch_index in 0..batches.batch_count() {
        let batch = batches.batch(batch_index);
        while cursors[batch_index] < batch.len() {
            let candidate = batch[cursors[batch_index]].clone();
            cursors[batch_index] += 1;
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

#[cfg(test)]
mod tests {
    use super::CandidateBatches;
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};

    fn candidate(id: &str) -> crate::domain::BlueprintCandidate {
        canonical_option2_candidate(id, operating_point(2.0e-6, 30_000.0, 0.18))
    }

    #[test]
    fn batches_keep_their_own_runs_and_cursors_start_at_zero() {
        let mut batches = CandidateBatches::new();
        batches.push_batch(vec![candidate("a"), candidate("b")]);
        batches.push_batch(Vec::new());
        batches.push_batch(vec![candidate("c")]);

        assert_eq!(batches.batch_count(), 3);
        assert_eq!(batches.batch(0).len(), 2);
        assert_eq!(batches.batch(0)[0].id, "a");
        assert_eq!(batches.batch(0)[1].id, "b");
        assert!(batches.batch(1).is_empty());
        assert_eq!(batches.batch(2)[0].id, "c");
        assert_eq!(batches.cursors(), vec![0, 0, 0]);
    }

    #[test]
    fn an_empty_set_has_no_batches_and_no_cursors() {
        let batches = CandidateBatches::new();
        assert_eq!(batches.batch_count(), 0);
        assert!(batches.cursors().is_empty());
    }

    #[test]
    fn results_assemble_in_order() {
        let batches = CandidateBatches::from_results(vec![
            Ok(vec![candidate("a")]),
            Ok(vec![candidate("b")]),
        ])
        .expect("invariant: fixture batches are Ok");
        assert_eq!(batches.batch_count(), 2);
        assert_eq!(batches.batch(1)[0].id, "b");
    }
}
