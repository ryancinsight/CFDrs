//! Phase 7 — Optimizer integration tests for `cfd-optim`.
//!
//! These tests exercise the full `SdtOptimizer` pipeline: candidate space
//! generation, physics evaluation, scoring, and ranking.

use cfd_optim::{OptimMode, SdtOptimizer, SdtWeights};

// ---------------------------------------------------------------------------
// Test 1: top_k(1) returns exactly one result
// ---------------------------------------------------------------------------

#[test]
fn top_k_1_returns_one_result() {
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let results = optimizer.top_k(1).unwrap();
    assert_eq!(results.len(), 1, "top_k(1) must return exactly 1 result");
}

// ---------------------------------------------------------------------------
// Test 2: top_k(5) returns exactly five results
// ---------------------------------------------------------------------------

#[test]
fn top_k_5_returns_five_results() {
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let results = optimizer.top_k(5).unwrap();
    assert_eq!(results.len(), 5, "top_k(5) must return exactly 5 results");
}

// ---------------------------------------------------------------------------
// Test 3: SdtCavitation top-1 has cavitation_number < 1 (active cavitation)
// ---------------------------------------------------------------------------

#[test]
fn sdt_cavitation_top1_has_sigma_below_one() {
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let results = optimizer.top_k(1).unwrap();
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
    let optimizer = SdtOptimizer::new(OptimMode::SdtTherapy, SdtWeights::default());
    let results = optimizer.top_k(5).unwrap();
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
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let results = optimizer.top_k(5).unwrap();
    for d in &results {
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
    let optimizer = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default());
    let results = optimizer.top_k(5).unwrap();
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
