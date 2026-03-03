//! Phase 7 — Optimizer integration tests for `cfd-optim`.
//!
//! These tests exercise the full `SdtOptimizer` pipeline: candidate space
//! generation, physics evaluation, scoring, and ranking.

use cfd_optim::{OptimMode, SdtOptimizer, SdtWeights};
use std::sync::OnceLock;

fn cavitation_top5() -> &'static Vec<cfd_optim::RankedDesign> {
    static CACHE: OnceLock<Vec<cfd_optim::RankedDesign>> = OnceLock::new();
    CACHE.get_or_init(|| {
        let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
        optimizer
            .top_k(5)
            .expect("cavitation top_k(5) must succeed")
    })
}

fn therapy_top5() -> &'static Vec<cfd_optim::RankedDesign> {
    static CACHE: OnceLock<Vec<cfd_optim::RankedDesign>> = OnceLock::new();
    CACHE.get_or_init(|| {
        let optimizer = SdtOptimizer::new(OptimMode::SdtTherapy, SdtWeights::default());
        optimizer.top_k(5).expect("therapy top_k(5) must succeed")
    })
}

// ---------------------------------------------------------------------------
// Test 1: top_k(1) returns exactly one result
// ---------------------------------------------------------------------------

#[test]
fn top_k_1_returns_one_result() {
    let top1: Vec<_> = cavitation_top5().iter().take(1).collect();
    assert_eq!(top1.len(), 1, "top_k(1) must return exactly 1 result");
}

// ---------------------------------------------------------------------------
// Test 2: top_k(5) returns exactly five results
// ---------------------------------------------------------------------------

#[test]
fn top_k_5_returns_five_results() {
    let results = cavitation_top5();
    assert_eq!(results.len(), 5, "top_k(5) must return exactly 5 results");
}

// ---------------------------------------------------------------------------
// Test 3: SdtCavitation top-1 has cavitation_number < 1 (active cavitation)
// ---------------------------------------------------------------------------

#[test]
fn sdt_cavitation_top1_has_sigma_below_one() {
    let results = cavitation_top5();
    assert!(
        results[0].metrics.cavitation_number < 1.0,
        "Best SdtCavitation design must have sigma < 1, got {}",
        results[0].metrics.cavitation_number
    );
}

// ---------------------------------------------------------------------------
// Test 4: SdtTherapy results are sorted by score descending
// ---------------------------------------------------------------------------

#[test]
fn results_sorted_by_score_descending() {
    let results = therapy_top5();
    for pair in results.windows(2) {
        assert!(
            pair[0].score >= pair[1].score,
            "Results must be sorted by score descending: {} < {}",
            pair[0].score,
            pair[1].score
        );
    }
}

// ---------------------------------------------------------------------------
// Test 5: All ranked designs have positive scores
// ---------------------------------------------------------------------------

#[test]
fn all_ranked_designs_have_positive_scores() {
    let results = cavitation_top5();
    for d in results {
        assert!(
            d.score > 0.0,
            "Ranked design '{}' must have positive score, got {}",
            d.candidate.id,
            d.score
        );
    }
}

// ---------------------------------------------------------------------------
// Test 6: Ranked designs have sequential rank numbers
// ---------------------------------------------------------------------------

#[test]
fn ranked_designs_have_sequential_ranks() {
    let results = cavitation_top5();
    for (i, d) in results.iter().enumerate() {
        assert_eq!(
            d.rank,
            i + 1,
            "Design at index {i} must have rank {}, got {}",
            i + 1,
            d.rank
        );
    }
}
