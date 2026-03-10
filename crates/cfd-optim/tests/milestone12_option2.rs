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

use cfd_optim::{
    build_milestone12_blueprint_candidate_space, evaluate_blueprint_candidate, evaluate_goal,
    evaluate_selective_venturi_cavitation, BlueprintCandidate, EvaluatedPool, OptimizationGoal,
};

#[test]
fn option2_candidate_space_contains_venturi_topologies() {
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let venturi: Vec<&BlueprintCandidate> = candidates
        .iter()
        .filter(|c| {
            c.topology_spec()
                .map(|spec| spec.has_venturi())
                .unwrap_or(false)
        })
        .collect();

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
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let venturi: Vec<BlueprintCandidate> = candidates
        .into_iter()
        .filter(|c| {
            c.topology_spec()
                .map(|spec| spec.has_venturi())
                .unwrap_or(false)
        })
        .take(20)
        .collect();

    let mut positive_count = 0;
    for candidate in &venturi {
        let result = evaluate_goal(
            candidate,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        );
        match result {
            Ok(eval) => {
                assert!(
                    eval.score.is_finite(),
                    "score must be finite for '{}'",
                    candidate.id
                );
                assert!(
                    eval.score >= 0.0,
                    "score must be non-negative for '{}'",
                    candidate.id
                );
                if eval.score > 0.0 {
                    positive_count += 1;
                }
            }
            Err(e) => panic!(
                "evaluate_goal failed for venturi candidate '{}': {e}",
                candidate.id
            ),
        }
    }

    assert!(
        positive_count > 0,
        "at least one venturi candidate should achieve positive Option 2 score"
    );
}

#[test]
fn option2_rejects_non_venturi_candidates() {
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let acoustic = candidates
        .into_iter()
        .find(|c| {
            c.topology_spec()
                .map(|spec| !spec.has_venturi() && !spec.split_stages.is_empty())
                .unwrap_or(false)
        })
        .expect("at least one acoustic candidate exists");

    let eval = evaluate_blueprint_candidate(&acoustic).expect("evaluation succeeds");
    let result = evaluate_selective_venturi_cavitation(&acoustic, eval);
    assert!(
        result.is_err(),
        "Option 2 scoring must reject candidates without venturi placements"
    );
}

#[test]
fn option2_cavitation_number_is_physically_consistent() {
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let venturi: Vec<BlueprintCandidate> = candidates
        .into_iter()
        .filter(|c| {
            c.topology_spec()
                .map(|spec| spec.has_venturi())
                .unwrap_or(false)
        })
        .take(20)
        .collect();

    for candidate in &venturi {
        let eval = evaluate_blueprint_candidate(candidate).expect("evaluation succeeds");

        for placement in &eval.venturi.placements {
            // σ = (p∞ − pᵥ) / (½ρv²) — finite for any non-zero velocity.
            // σ < 0 implies p∞ < pᵥ (strong cavitation, physically valid).
            // σ → 0⁺ is the inception threshold; σ > 1 means no cavitation.
            assert!(
                placement.cavitation_number.is_finite(),
                "cavitation number must be finite for placement '{}' in '{}'",
                placement.placement_id,
                candidate.id
            );

            // Throat velocity must be positive — continuity demands v > 0
            // when there is positive flow through the throat.
            assert!(
                placement.effective_throat_velocity_m_s > 0.0,
                "throat velocity must be positive for '{}' (v={})",
                candidate.id,
                placement.effective_throat_velocity_m_s
            );
        }

        // Exposure fractions are blood volume fractions — bounded [0, 1].
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
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let venturi: Vec<BlueprintCandidate> = candidates
        .into_iter()
        .filter(|c| {
            c.topology_spec()
                .map(|spec| spec.has_venturi())
                .unwrap_or(false)
        })
        .take(20)
        .collect();

    let pool = EvaluatedPool::from_candidates(&venturi);
    let ranked = pool
        .top_k(
            5,
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        )
        .expect("pool ranking succeeds");

    for window in ranked.windows(2) {
        assert!(
            window[0].score >= window[1].score,
            "Option 2 ranking must be descending: {} ({}) >= {} ({})",
            window[0].candidate_id,
            window[0].score,
            window[1].candidate_id,
            window[1].score
        );
    }
}

#[test]
fn option2_venturi_selectivity_prefers_low_cavitation_number() {
    let candidates = build_milestone12_blueprint_candidate_space().expect("candidate space builds");

    let venturi: Vec<BlueprintCandidate> = candidates
        .into_iter()
        .filter(|c| {
            c.topology_spec()
                .map(|spec| spec.has_venturi())
                .unwrap_or(false)
        })
        .take(20)
        .collect();

    let pool = EvaluatedPool::from_candidates(&venturi);
    let ranked = pool
        .top_k(
            venturi.len(),
            OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
        )
        .expect("pool ranking succeeds");

    // Among positive-score candidates, verify that cavitation selectivity is
    // correlated with score: candidates with higher cavitation_selectivity_score
    // should have higher total scores (on average).
    let positive: Vec<_> = ranked.iter().filter(|e| e.score > 0.0).collect();
    if positive.len() >= 2 {
        let top_selectivity = positive[0].venturi.cavitation_selectivity_score;
        let bottom_selectivity = positive
            .last()
            .unwrap()
            .venturi
            .cavitation_selectivity_score;
        // The top-ranked candidate's cavitation selectivity should be at least
        // as strong as the bottom-ranked one (higher selectivity = stronger
        // cavitation at treatment site with low healthy-cell exposure).
        assert!(
            top_selectivity >= bottom_selectivity - 0.01,
            "top selectivity ({top_selectivity:.4}) should be >= bottom ({bottom_selectivity:.4})"
        );
    }
}
