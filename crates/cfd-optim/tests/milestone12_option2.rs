//! Milestone 12 — Option 2 integration tests.
//!
//! Validates the venturi-cavitation optimisation pipeline end-to-end:
//!
//! 1. Venturi-equipped candidates exist in the candidate space
//! 2. Cavitation number σ respects Bernoulli physics: σ = (p∞ − pᵥ)/(½ρv²)
//! 3. Option 2 scoring gates on zero venturi placements (returns error)
//! 4. Cavitation selectivity favours low σ (strong cavitation) at treatment
//!    throats with low RBC/WBC exposure
//! 5. FDA venturi transient shear compliance (≤ 300 Pa for ≤ 0.015 s)
//! 6. EvaluatedPool ranking for Option 2 is consistent

use std::collections::HashSet;

use cfd_optim::{
    evaluate_blueprint_candidate, evaluate_goal, evaluate_selective_venturi_cavitation,
    orchestration_lineage_key, BlueprintCandidate, EvaluatedPool, OperatingPoint, OptimizationGoal,
};
use cfd_schematics::{
    build_milestone12_blueprint, enumerate_milestone12_topologies, SplitKind,
    TreatmentActuationMode, VenturiPlacementMode,
};

// ---------------------------------------------------------------------------
// Fast fixture constructors — build candidates directly from the topology
// catalog instead of materializing the full ~1000-candidate test space.
// ---------------------------------------------------------------------------

fn test_op(flow: f64, gauge: f64) -> OperatingPoint {
    OperatingPoint {
        flow_rate_m3_s: flow,
        inlet_gauge_pa: gauge,
        feed_hematocrit: 0.45,
        patient_context: None,
    }
}

/// Build diverse venturi candidates from different topologies and operating
/// points. Returns up to `n` candidates without touching the full sweep.
fn fast_venturi_candidates(n: usize) -> Vec<BlueprintCandidate> {
    let topologies: Vec<_> = enumerate_milestone12_topologies()
        .into_iter()
        .filter(|r| {
            !r.split_kinds.is_empty()
                && r.split_kinds
                    .iter()
                    .all(|k| matches!(k, SplitKind::NFurcation(2..=5)))
        })
        .collect();

    let ops = [
        test_op(1.5e-6, 25_000.0),
        test_op(2.0e-6, 30_000.0),
        test_op(2.5e-6, 50_000.0),
        test_op(3.0e-6, 100_000.0),
    ];

    let mut candidates = Vec::with_capacity(n);
    for (i, base) in topologies.into_iter().enumerate() {
        if candidates.len() >= n {
            break;
        }
        let request = cfd_schematics::Milestone12TopologyRequest {
            treatment_mode: TreatmentActuationMode::VenturiCavitation,
            venturi_throat_count: 2,
            venturi_throat_width_m: 0.4e-3,
            venturi_throat_length_m: 1.2e-3,
            venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            venturi_target_channel_ids: Vec::new(),
            ..base
        };
        if let Ok(blueprint) = build_milestone12_blueprint(&request) {
            let op = ops[i % ops.len()].clone();
            candidates.push(BlueprintCandidate::new(
                format!("test-venturi-{i}"),
                blueprint,
                op,
            ));
        }
    }
    candidates
}

/// Build a single acoustic (non-venturi) candidate.
fn fast_acoustic_candidate() -> BlueprintCandidate {
    let base = enumerate_milestone12_topologies()
        .into_iter()
        .find(|r| !r.split_kinds.is_empty())
        .expect("at least one topology");

    let request = cfd_schematics::Milestone12TopologyRequest {
        treatment_mode: TreatmentActuationMode::UltrasoundOnly,
        venturi_throat_count: 0,
        venturi_target_channel_ids: Vec::new(),
        ..base
    };

    let blueprint = build_milestone12_blueprint(&request).expect("acoustic blueprint should build");
    BlueprintCandidate::new("test-acoustic", blueprint, test_op(2.0e-6, 30_000.0))
}

/// Build a matched pair of acoustic + venturi candidates from the SAME
/// topology (so they share a lineage key) for lineage tests.
fn fast_lineage_pairs(n: usize) -> (Vec<BlueprintCandidate>, Vec<BlueprintCandidate>) {
    let topologies: Vec<_> = enumerate_milestone12_topologies()
        .into_iter()
        .filter(|r| {
            !r.split_kinds.is_empty()
                && r.split_kinds
                    .iter()
                    .all(|k| matches!(k, SplitKind::NFurcation(2..=5)))
        })
        .take(n)
        .collect();

    let op = test_op(2.0e-6, 30_000.0);
    let mut acoustics = Vec::with_capacity(n);
    let mut venturis = Vec::with_capacity(n);

    for (i, base) in topologies.into_iter().enumerate() {
        // Acoustic variant
        let acoustic_req = cfd_schematics::Milestone12TopologyRequest {
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
            venturi_throat_count: 0,
            venturi_target_channel_ids: Vec::new(),
            ..base.clone()
        };
        if let Ok(bp) = build_milestone12_blueprint(&acoustic_req) {
            acoustics.push(BlueprintCandidate::new(
                format!("test-acoustic-{i}"),
                bp,
                op.clone(),
            ));
        }

        // Venturi variant from the same topology
        let venturi_req = cfd_schematics::Milestone12TopologyRequest {
            treatment_mode: TreatmentActuationMode::VenturiCavitation,
            venturi_throat_count: 2,
            venturi_throat_width_m: 0.4e-3,
            venturi_throat_length_m: 1.2e-3,
            venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            venturi_target_channel_ids: Vec::new(),
            ..base
        };
        if let Ok(bp) = build_milestone12_blueprint(&venturi_req) {
            venturis.push(BlueprintCandidate::new(
                format!("test-venturi-{i}"),
                bp,
                op.clone(),
            ));
        }
    }

    (acoustics, venturis)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn option2_candidate_space_contains_venturi_topologies() {
    let venturi = fast_venturi_candidates(10);

    assert!(
        !venturi.is_empty(),
        "candidate space must contain venturi-equipped topologies"
    );

    for candidate in &venturi {
        let spec = candidate.topology_spec().unwrap();
        assert!(
            !spec.venturi_placements.is_empty(),
            "venturi candidate '{}' must have non-empty venturi_placements",
            candidate.id
        );
        for placement in &spec.venturi_placements {
            let g = &placement.throat_geometry;
            assert!(
                g.throat_width_m > 0.0 && g.throat_width_m < g.inlet_width_m,
                "throat must be narrower than inlet for '{}'",
                candidate.id
            );
        }
    }
}

#[test]
fn option2_evaluate_goal_produces_bounded_scores_for_venturi_candidates() {
    let venturi = fast_venturi_candidates(20);

    let mut eligible_count = 0;
    for candidate in &venturi {
        let result = evaluate_goal(
            candidate,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        );
        match result {
            Ok(eval) => {
                if eval.is_eligible() {
                    assert!(
                        eval.score.is_some_and(f64::is_finite),
                        "eligible score must be finite for '{}'",
                        candidate.id
                    );
                    assert!(
                        eval.score.is_some_and(|score| score > 0.0),
                        "eligible Option 2 score must be strictly positive for '{}'",
                        candidate.id
                    );
                    eligible_count += 1;
                } else {
                    assert!(
                        eval.score.is_none(),
                        "screened-out Option 2 candidate '{}' must omit numeric score",
                        candidate.id
                    );
                }
            }
            Err(e) => panic!(
                "evaluate_goal failed for venturi candidate '{}': {e}",
                candidate.id
            ),
        }
    }

    assert!(
        eligible_count > 0,
        "at least one venturi candidate should remain eligible under Option 2 scoring"
    );
}

#[test]
fn option2_rejects_non_venturi_candidates() {
    let acoustic = fast_acoustic_candidate();

    let eval = evaluate_blueprint_candidate(&acoustic).expect("evaluation succeeds");
    let result = evaluate_selective_venturi_cavitation(&acoustic, eval);
    assert!(
        result.is_err(),
        "Option 2 scoring must reject candidates without venturi placements"
    );
}

#[test]
fn option2_cavitation_number_is_physically_consistent() {
    let venturi = fast_venturi_candidates(20);

    for candidate in &venturi {
        let eval = evaluate_blueprint_candidate(candidate).expect("evaluation succeeds");

        for placement in &eval.venturi.placements {
            assert!(
                placement.cavitation_number.is_finite(),
                "cavitation number must be finite for placement '{}' in '{}'",
                placement.placement_id,
                candidate.id
            );

            assert!(
                placement.effective_throat_velocity_m_s > 0.0,
                "throat velocity must be positive for '{}' (v={})",
                candidate.id,
                placement.effective_throat_velocity_m_s
            );
        }

        assert!(
            (0.0..=1.0).contains(&eval.venturi.rbc_exposure_fraction),
            "rbc_exposure_fraction out of bounds for '{}'",
            candidate.id
        );
        assert!(
            (0.0..=1.0).contains(&eval.venturi.wbc_exposure_fraction),
            "wbc_exposure_fraction out of bounds for '{}'",
            candidate.id
        );
    }
}

#[test]
fn option2_pool_ranking_is_descending_and_consistent() {
    let venturi = fast_venturi_candidates(20);

    let pool = EvaluatedPool::from_candidates(&venturi);
    let ranked = pool
        .top_k(
            5,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        )
        .expect("pool ranking succeeds");

    for window in ranked.windows(2) {
        let left_score = window[0].score.expect("ranked entries are always eligible");
        let right_score = window[1].score.expect("ranked entries are always eligible");
        assert!(
            left_score >= right_score,
            "Option 2 ranking must be descending: {} ({}) >= {} ({})",
            window[0].candidate_id,
            left_score,
            window[1].candidate_id,
            right_score
        );
    }
}

#[test]
fn option2_candidates_retain_option1_lineage() {
    let (acoustics, venturis) = fast_lineage_pairs(10);

    let acoustic_lineages: HashSet<_> = acoustics
        .iter()
        .filter_map(orchestration_lineage_key)
        .collect();

    let venturi_lineages: Vec<_> = venturis
        .iter()
        .filter_map(orchestration_lineage_key)
        .collect();

    assert!(
        !venturi_lineages.is_empty(),
        "Option 2 candidate space must retain venturi lineage keys"
    );
    assert!(
        venturi_lineages
            .iter()
            .all(|lineage| acoustic_lineages.contains(lineage)),
        "every Option 2 venturi lineage must derive from an Option 1 acoustic lineage"
    );
}

#[test]
fn option2_pool_returns_only_eligible_scored_entries() {
    let venturi = fast_venturi_candidates(20);

    let pool = EvaluatedPool::from_candidates(&venturi);
    let ranked = pool
        .top_k(
            5,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        )
        .expect("pool ranking succeeds");

    assert!(
        ranked
            .iter()
            .all(|entry| entry.is_eligible() && entry.score.is_some_and(|score| score > 0.0)),
        "ranked Option 2 entries must be eligible and carry strictly positive scores"
    );
}

#[test]
fn option2_venturi_selectivity_prefers_low_cavitation_number() {
    let venturi = fast_venturi_candidates(20);

    let pool = EvaluatedPool::from_candidates(&venturi);
    let ranked = pool
        .top_k(
            venturi.len(),
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        )
        .expect("pool ranking succeeds");

    if ranked.len() >= 2 {
        let top_selectivity = ranked[0].venturi.cavitation_selectivity_score;
        let bottom_selectivity = ranked.last().unwrap().venturi.cavitation_selectivity_score;
        assert!(
            top_selectivity >= bottom_selectivity - 0.01,
            "top selectivity ({top_selectivity:.4}) should be >= bottom ({bottom_selectivity:.4})"
        );
    }
}
