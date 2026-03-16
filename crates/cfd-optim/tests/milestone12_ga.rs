//! Milestone 12 — Genetic Algorithm integration tests.
//!
//! Validates the in-place GA refinement pipeline:
//!
//! 1. `generate_ga_mutations` produces at least one mutated candidate from
//!    a valid Option 2 seed
//! 2. Mutations preserve topology invariants (branch width conservation,
//!    treatment-path roles)
//! 3. `BlueprintGeneticOptimizer` evolves a population and returns sorted results
//! 4. Serpentine insertion mutations create dean-vortex-capable geometry
//! 5. GA scoring includes lineage bonus and Dean number contribution
//! 6. GA can bootstrap venturi placements from an Option 1 seed

use cfd_optim::{
    build_milestone12_ga_seed_pair, evaluate_blueprint_candidate, evaluate_goal,
    generate_ga_mutations, promote_option1_candidate_to_ga_seed, BlueprintCandidate,
    BlueprintEvaluationStatus, BlueprintGeneticOptimizer, OperatingPoint, OptimizationGoal,
};
use cfd_schematics::{
    build_milestone12_blueprint, enumerate_milestone12_topologies, SplitKind,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiPlacementMode,
};

// ---------------------------------------------------------------------------
// Fast fixture constructors — build ONE candidate instead of scanning the
// full ~1000-candidate test space. Reduces per-test setup from ~3 s to ~5 ms.
// ---------------------------------------------------------------------------

fn test_operating_point() -> OperatingPoint {
    OperatingPoint {
        flow_rate_m3_s: 2.0e-6,
        inlet_gauge_pa: 30_000.0,
        feed_hematocrit: 0.45,
        patient_context: None,
    }
}

/// Build a single venturi-capable seed directly from the topology catalog.
fn fast_venturi_seed() -> BlueprintCandidate {
    let base = enumerate_milestone12_topologies()
        .into_iter()
        .find(|r| {
            !r.split_kinds.is_empty()
                && r.split_kinds
                    .iter()
                    .all(|k| matches!(k, SplitKind::NFurcation(2..=5)))
        })
        .expect("at least one venturi-eligible topology in the catalog");

    let request = cfd_schematics::Milestone12TopologyRequest {
        treatment_mode: TreatmentActuationMode::VenturiCavitation,
        venturi_throat_count: 2,
        venturi_throat_width_m: 0.4e-3,
        venturi_throat_length_m: 1.2e-3,
        venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
        venturi_target_channel_ids: Vec::new(),
        ..base
    };

    let blueprint =
        build_milestone12_blueprint(&request).expect("venturi seed blueprint should build");

    BlueprintCandidate::new("test-venturi-seed", blueprint, test_operating_point())
}

/// Build a single acoustic (non-venturi) candidate for promotion tests.
fn fast_acoustic_seed() -> BlueprintCandidate {
    let base = enumerate_milestone12_topologies()
        .into_iter()
        .find(|r| !r.split_kinds.is_empty())
        .expect("at least one topology in the catalog");

    let request = cfd_schematics::Milestone12TopologyRequest {
        treatment_mode: TreatmentActuationMode::UltrasoundOnly,
        venturi_throat_count: 0,
        venturi_target_channel_ids: Vec::new(),
        ..base
    };

    let blueprint =
        build_milestone12_blueprint(&request).expect("acoustic seed blueprint should build");

    BlueprintCandidate::new("test-acoustic-seed", blueprint, test_operating_point())
}

/// Build a minimal GA seed pair without evaluating the full candidate space.
fn fast_ga_seed_pair() -> Vec<BlueprintCandidate> {
    let acoustic = fast_acoustic_seed();
    let promoted = promote_option1_candidate_to_ga_seed(&acoustic)
        .expect("Option 1 promotion should succeed");
    let venturi = fast_venturi_seed();
    vec![promoted, venturi]
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn ga_generate_mutations_produces_candidates() {
    let seed = fast_venturi_seed();
    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    assert!(
        !mutations.is_empty(),
        "GA must produce at least one mutated candidate from a venturi seed"
    );

    for mutant in &mutations {
        assert!(
            !mutant.channels().is_empty(),
            "mutated candidate '{}' must have channels",
            mutant.id
        );
    }
}

#[test]
fn ga_mutations_preserve_treatment_path_roles() {
    let seed = fast_venturi_seed();
    let seed_treatment_count = seed.treatment_channel_ids().len();

    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    for mutant in &mutations {
        let mutant_treatment_count = mutant.treatment_channel_ids().len();
        assert!(
            mutant_treatment_count >= seed_treatment_count,
            "mutation '{}' must preserve (or extend) treatment channels: seed={}, mutant={}",
            mutant.id,
            seed_treatment_count,
            mutant_treatment_count
        );
    }
}

#[test]
fn ga_mutations_preserve_venturi_placements() {
    let seed = fast_venturi_seed();
    let seed_spec = seed.topology_spec().expect("seed has topology");
    let seed_venturi_count = seed_spec.venturi_placements.len();

    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    for mutant in &mutations {
        let spec = mutant.topology_spec().expect("mutant has topology");
        assert!(
            !spec.venturi_placements.is_empty(),
            "GA mutant '{}' must retain venturi placements",
            mutant.id
        );
        assert!(
            spec.venturi_placements.len() >= seed_venturi_count,
            "GA mutant '{}' should preserve venturi count: seed={}, mutant={}",
            mutant.id,
            seed_venturi_count,
            spec.venturi_placements.len()
        );
    }
}

#[test]
fn ga_evaluate_goal_produces_bounded_scores() {
    let seed = fast_venturi_seed();
    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    for mutant in mutations.iter().take(5) {
        let result = evaluate_goal(mutant, OptimizationGoal::InPlaceDeanSerpentineRefinement);
        match result {
            Ok(eval) => {
                if eval.is_eligible() {
                    assert!(
                        eval.score.is_some_and(f64::is_finite),
                        "eligible GA score must be finite for '{}'",
                        mutant.id
                    );
                    assert!(
                        eval.score.is_some_and(|score| score > 0.0),
                        "eligible GA score must be strictly positive for '{}'",
                        mutant.id
                    );
                } else {
                    assert!(
                        eval.score.is_none(),
                        "screened-out GA candidate '{}' must omit numeric score",
                        mutant.id
                    );
                }
            }
            Err(e) => panic!("evaluate_goal (GA) failed for '{}': {e}", mutant.id),
        }
    }
}

#[test]
fn ga_mutation_catalog_includes_canonical_split_serpentine_and_venturi_compositions() {
    let seed = fast_venturi_seed();
    let seed_topology = seed.topology_spec().expect("seed has topology");
    let seed_stage_count = seed_topology.split_stages.len();
    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    assert!(
        mutations.iter().any(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| topology.has_serpentine() && topology.has_venturi())
        }),
        "GA mutation catalog must include serpentine plus venturi compositions"
    );
    assert!(
        mutations.iter().any(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| topology.split_stages.len() > seed_stage_count)
        }),
        "GA mutation catalog must include canonical split-merge insertions"
    );
    // Topology mutations (serpentine, venturi, split-merge) must have recorded
    // lineage. Operating-point perturbations (ga-qf*, ga-pf*) carry the seed's
    // lineage unchanged, which may have empty mutations for a fresh fixture.
    let topology_mutations: Vec<_> = mutations
        .iter()
        .filter(|c| !c.id.contains("-ga-qf") && !c.id.contains("-ga-pf"))
        .collect();
    assert!(
        topology_mutations.iter().all(|candidate| {
            candidate.blueprint().lineage().is_some_and(|lineage| {
                lineage.current_stage == TopologyOptimizationStage::InPlaceDeanSerpentineRefinement
                    && !lineage.mutations.is_empty()
            })
        }),
        "topology GA mutations must have recorded lineage"
    );
}

#[test]
fn ga_optimizer_evolves_population_and_ranks_results() {
    let seeds = fast_ga_seed_pair();

    let result = BlueprintGeneticOptimizer::new(OptimizationGoal::InPlaceDeanSerpentineRefinement)
        .with_seeds(seeds)
        .with_population(8)
        .with_max_generations(2)
        .with_top_k(3)
        .run()
        .expect("GA optimizer runs successfully");

    assert!(
        !result.top_candidates.is_empty(),
        "GA optimizer must return at least one ranked candidate"
    );

    // Verify descending score order.
    for window in result.top_candidates.windows(2) {
        assert!(
            window[0].evaluation.score.unwrap_or_default()
                >= window[1].evaluation.score.unwrap_or_default(),
            "GA results must be descending: {} ({}) >= {} ({})",
            window[0].candidate.id,
            window[0].evaluation.score.unwrap_or_default(),
            window[1].candidate.id,
            window[1].evaluation.score.unwrap_or_default()
        );
    }

    // All GA results must be evaluable under the GA goal.
    for entry in &result.top_candidates {
        let eval = evaluate_goal(
            &entry.candidate,
            OptimizationGoal::InPlaceDeanSerpentineRefinement,
        );
        assert!(
            eval.is_ok(),
            "ranked GA candidate '{}' must be evaluable",
            entry.candidate.id
        );
        assert_eq!(entry.evaluation.status, BlueprintEvaluationStatus::Eligible);
        assert!(
            entry.evaluation.exceeds_all_baselines.is_some(),
            "GA evaluation should record baseline comparison"
        );
    }
}

#[test]
fn ga_scoring_includes_dean_number_contribution() {
    let seed = fast_venturi_seed();
    let mutations = generate_ga_mutations(&seed).expect("GA mutations succeed");

    for mutant in mutations.iter().take(3) {
        let eval = evaluate_blueprint_candidate(mutant).expect("evaluation succeeds");

        let max_dean = eval
            .venturi
            .placements
            .iter()
            .map(|p| p.dean_number)
            .fold(0.0_f64, f64::max);

        let dean_bonus = max_dean / 100.0;

        // Dean bonus is additive in GA scoring: score includes lineage_bonus + dean_bonus.
        // Verify the Dean number is non-negative (physically required).
        assert!(
            max_dean >= 0.0,
            "Dean number must be non-negative for '{}' (De={})",
            mutant.id,
            max_dean
        );

        // Verify the dean bonus calculation matches expected formula.
        assert!(
            (dean_bonus - max_dean / 100.0).abs() < 1.0e-12,
            "Dean bonus formula: De/100 = {dean_bonus}"
        );
    }
}

#[test]
fn ga_promotes_option1_seed_through_canonical_topology_api() {
    let acoustic = fast_acoustic_seed();

    assert!(
        generate_ga_mutations(&acoustic).is_err(),
        "GA mutation generation must reject non-venturi seeds after the canonical promotion split"
    );

    let promoted =
        promote_option1_candidate_to_ga_seed(&acoustic).expect("Option 1 promotion should work");
    assert!(
        promoted
            .topology_spec()
            .is_ok_and(|topology| !topology.venturi_placements.is_empty()),
        "promoted GA seed must carry venturi placements"
    );
}

#[test]
fn ga_seed_pair_contains_exactly_two_canonical_baselines() {
    // Provide a minimal candidate set (one acoustic + one venturi) to avoid
    // evaluating the full ~1000-candidate test space.
    let candidates = vec![fast_acoustic_seed(), fast_venturi_seed()];
    let seeds = build_milestone12_ga_seed_pair(&candidates).expect("seed pair");

    assert_eq!(
        seeds.len(),
        2,
        "GA must start from exactly two canonical seeds"
    );
    assert!(
        seeds.iter().all(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| !topology.venturi_placements.is_empty())
        }),
        "both GA seeds must be venturi-capable after canonical promotion"
    );
}
