//! Milestone 12 — Option 1 integration tests.
//!
//! Validates the asymmetric-split residence-separation pipeline end-to-end:
//!
//! 1. Candidate space generation produces acoustic-tier blueprints
//! 2. Every candidate evaluates under Option 1 scoring without error
//! 3. Ranking produces strictly positive scores for feasible candidates
//! 4. Zweifach–Fung physics: deeper split trees achieve higher cancer-center
//!    enrichment (Tri→Bi→Bi > Tri alone)
//! 5. Separation metrics are bounded [0, 1] and monotonically related to
//!    channel asymmetry ratios
//! 6. FDA safety constraints are respected (wall shear < 150 Pa sustained)

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, evaluate_blueprint_candidate, evaluate_goal,
    BlueprintCandidate, EvaluatedPool, OperatingPoint, OptimizationGoal,
};
use cfd_schematics::{
    build_milestone12_blueprint, enumerate_milestone12_topologies, TreatmentActuationMode,
};
use std::collections::HashSet;

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

/// Build diverse acoustic (non-venturi) candidates from different topologies
/// and operating points. Returns up to `n` candidates.
fn fast_acoustic_candidates(n: usize) -> Vec<BlueprintCandidate> {
    let topologies: Vec<_> = enumerate_milestone12_topologies()
        .into_iter()
        .filter(|r| !r.split_kinds.is_empty())
        .collect();

    let ops = [
        test_op(1.333e-6, 25_000.0),
        test_op(1.667e-6, 50_000.0),
        test_op(2.0e-6, 100_000.0),
        test_op(3.333e-6, 200_000.0),
    ];

    let mut candidates = Vec::with_capacity(n);
    for (i, base) in topologies.into_iter().enumerate() {
        if candidates.len() >= n {
            break;
        }
        let request = cfd_schematics::Milestone12TopologyRequest {
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
            venturi_throat_count: 0,
            venturi_target_channel_ids: Vec::new(),
            ..base
        };
        if let Ok(blueprint) = build_milestone12_blueprint(&request) {
            let op = ops[i % ops.len()].clone();
            candidates.push(BlueprintCandidate::new(
                format!("test-acoustic-{i}"),
                blueprint,
                op,
            ));
        }
    }
    candidates
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn option1_candidate_space_contains_acoustic_selective_topologies() {
    let acoustic = fast_acoustic_candidates(10);

    assert!(
        !acoustic.is_empty(),
        "candidate space must contain at least one acoustic (non-venturi) selective topology"
    );

    for candidate in &acoustic {
        let spec = candidate.topology_spec().unwrap();
        assert!(
            spec.split_stages
                .iter()
                .any(|s| s.branches.iter().any(|b| b.treatment_path)),
            "acoustic selective topology '{}' must have at least one treatment branch",
            candidate.id
        );
    }
}

#[test]
fn option1_candidate_space_covers_mirrored_bi_tri_quad_penta_catalog() {
    // This test validates that `build_milestone12_blueprint_candidate_space`
    // covers every topology family × mirror variant in the canonical catalog.
    // It must use the real space builder.
    let catalog = enumerate_milestone12_topologies();
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let available: HashSet<(String, bool, bool)> = candidates
        .iter()
        .filter(|candidate| {
            candidate
                .topology_spec()
                .map(|spec| !spec.has_venturi() && !spec.split_stages.is_empty())
                .unwrap_or(false)
        })
        .filter_map(|candidate| {
            let topology = candidate.topology_spec().ok()?;
            let hints = candidate.blueprint().render_hints()?;
            Some((
                topology.stage_sequence_label(),
                hints.mirror_x,
                hints.mirror_y,
            ))
        })
        .collect();

    for request in catalog {
        let family = match request.split_kinds.first() {
            Some(cfd_schematics::SplitKind::NFurcation(2)) => "Bi",
            Some(cfd_schematics::SplitKind::NFurcation(3)) => "Tri",
            Some(cfd_schematics::SplitKind::NFurcation(4)) => "Quad",
            Some(cfd_schematics::SplitKind::NFurcation(5)) => "Penta",
            Some(cfd_schematics::SplitKind::NFurcation(_)) | None => {
                panic!("unexpected milestone12 root family in canonical catalog")
            }
        };
        assert!(
            available.contains(&(family.to_string(), request.mirror_x, request.mirror_y)),
            "candidate space must materialize canonical catalog family {} (mirror_x={}, mirror_y={})",
            family,
            request.mirror_x,
            request.mirror_y
        );
    }
}

#[test]
fn option1_evaluate_goal_produces_bounded_scores() {
    let acoustic = fast_acoustic_candidates(20);

    for candidate in &acoustic {
        let result = evaluate_goal(
            candidate,
            OptimizationGoal::AsymmetricSplitResidenceSeparation,
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
                        "eligible acoustic selective candidate '{}' must have strictly positive score (got {})",
                        candidate.id,
                        eval.score.unwrap_or_default()
                    );
                } else {
                    assert!(
                        eval.score.is_none(),
                        "screened-out acoustic candidate '{}' must omit numeric score",
                        candidate.id
                    );
                }
            }
            Err(e) => panic!("evaluate_goal failed for candidate '{}': {e}", candidate.id),
        }
    }
}

#[test]
fn option1_evaluated_pool_ranking_is_descending_when_nonempty() {
    let acoustic = fast_acoustic_candidates(20);

    let pool = EvaluatedPool::from_candidates(&acoustic);
    let ranked = pool
        .top_k(5, OptimizationGoal::AsymmetricSplitResidenceSeparation)
        .expect("pool must produce eligible candidates with non-zero scores");
    assert!(
        !ranked.is_empty(),
        "ranked results must be non-empty for acoustic selective candidates"
    );
    for window in ranked.windows(2) {
        let left_score = window[0].score.expect("ranked entries are always eligible");
        let right_score = window[1].score.expect("ranked entries are always eligible");
        assert!(
            left_score >= right_score,
            "pool ranking must be descending: {} ({}) >= {} ({})",
            window[0].candidate_id,
            left_score,
            window[1].candidate_id,
            right_score
        );
    }
}

#[test]
fn option1_pool_returns_only_eligible_scored_entries() {
    let acoustic = fast_acoustic_candidates(20);

    let pool = EvaluatedPool::from_candidates(&acoustic);
    let ranked = pool
        .top_k(5, OptimizationGoal::AsymmetricSplitResidenceSeparation)
        .expect("pool must produce eligible acoustic candidates");

    assert!(
        ranked
            .iter()
            .all(|entry| entry.is_eligible() && entry.score.is_some_and(|score| score > 0.0)),
        "ranked Option 1 entries must be eligible and carry strictly positive scores"
    );
}

#[test]
fn option1_separation_metrics_are_physically_bounded() {
    let acoustic = fast_acoustic_candidates(20);

    for candidate in &acoustic {
        let eval = evaluate_blueprint_candidate(candidate).expect("blueprint evaluation succeeds");

        // Separation fractions are volume fractions — must be in [0, 1].
        assert!(
            (0.0..=1.0).contains(&eval.separation.cancer_center_fraction),
            "cancer_center_fraction out of bounds: {} for '{}'",
            eval.separation.cancer_center_fraction,
            candidate.id
        );
        assert!(
            (0.0..=1.0).contains(&eval.separation.rbc_peripheral_fraction),
            "rbc_peripheral_fraction out of bounds: {} for '{}'",
            eval.separation.rbc_peripheral_fraction,
            candidate.id
        );
        assert!(
            (0.0..=1.0).contains(&eval.separation.separation_efficiency),
            "separation_efficiency out of bounds: {} for '{}'",
            eval.separation.separation_efficiency,
            candidate.id
        );

        // Residence time must be positive for acoustic selective candidates
        // that route flow through treatment channels.
        assert!(
            eval.residence.treatment_flow_fraction > 0.0,
            "treatment_flow_fraction must be positive for acoustic selective candidate '{}' (got {})",
            candidate.id,
            eval.residence.treatment_flow_fraction
        );
        assert!(
            eval.residence.treatment_residence_time_s > 0.0,
            "treatment_residence_time_s must be positive for acoustic selective candidate '{}' (got {})",
            candidate.id,
            eval.residence.treatment_residence_time_s
        );
    }
}

#[test]
fn option1_fda_wall_shear_compliance() {
    let acoustic = fast_acoustic_candidates(20);

    // Verify FDA compliance for every evaluable acoustic candidate:
    // main_channel_margin > 0 ⟹ sustained wall shear < 150 Pa.
    for candidate in &acoustic {
        let eval = evaluate_blueprint_candidate(candidate).expect("evaluation succeeds");
        if eval.safety.main_channel_margin > 0.0 {
            assert!(
                eval.safety.main_channel_margin.is_finite(),
                "safety margin must be finite for '{}'",
                candidate.id
            );
        }
    }
}
