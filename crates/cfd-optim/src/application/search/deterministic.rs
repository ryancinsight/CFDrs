use crate::application::objectives::{evaluate_goal, BlueprintObjectiveEvaluation};
use crate::domain::{BlueprintCandidate, OperatingPoint, OptimizationGoal};
use crate::error::OptimError;
use cfd_schematics::BlueprintTopologySpec;

pub fn build_blueprint_candidates_from_specs(
    specs: &[BlueprintTopologySpec],
    operating_point: &OperatingPoint,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    specs
        .iter()
        .enumerate()
        .map(|(index, spec)| {
            BlueprintCandidate::from_topology_spec(
                format!("{index:04}-{}", spec.topology_id),
                spec,
                operating_point.clone(),
            )
        })
        .collect()
}

pub fn rank_blueprint_candidates(
    goal: OptimizationGoal,
    candidates: &[BlueprintCandidate],
) -> Result<Vec<BlueprintObjectiveEvaluation>, OptimError> {
    let mut evaluations = candidates
        .iter()
        .map(|candidate| evaluate_goal(candidate, goal))
        .collect::<Result<Vec<_>, _>>()?;
    evaluations.sort_by(|left, right| {
        right
            .score
            .partial_cmp(&left.score)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| left.candidate_id.cmp(&right.candidate_id))
    });
    Ok(evaluations)
}

#[cfg(test)]
mod tests {
    use crate::domain::fixtures::{
        canonical_option1_spec, canonical_option2_candidate, operating_point,
    };

    use super::{build_blueprint_candidates_from_specs, rank_blueprint_candidates};

    #[test]
    fn official_option1_option2_ga_blueprints_are_evaluable_without_topology_enum() {
        let operating_point = operating_point(2.4e-6, 32_000.0, 0.12);
        let option1_candidates =
            build_blueprint_candidates_from_specs(&[canonical_option1_spec()], &operating_point)
                .expect("option1 candidates");
        let option2 = canonical_option2_candidate("option2", operating_point.clone());
        let ga = crate::generate_ga_mutations(&option2)
            .expect("ga mutations")
            .into_iter()
            .next()
            .expect("at least one ga candidate");

        let option1_ranked = rank_blueprint_candidates(
            crate::OptimizationGoal::SelectiveAcousticResidenceSeparation,
            &option1_candidates,
        )
        .expect("option1 evaluation");
        let option2_ranked = rank_blueprint_candidates(
            crate::OptimizationGoal::SelectiveVenturiCavitation,
            &[option2.clone()],
        )
        .expect("option2 evaluation");
        let ga_ranked =
            rank_blueprint_candidates(crate::OptimizationGoal::BlueprintGeneticRefinement, &[ga])
                .expect("ga evaluation");

        assert_eq!(option1_ranked.len(), 1);
        assert_eq!(option2_ranked.len(), 1);
        assert_eq!(ga_ranked.len(), 1);
    }

    #[test]
    fn canonical_series_specs_build_blueprint_candidates_without_bridge_logic() {
        let operating_point = operating_point(2.4e-6, 32_000.0, 0.12);
        let specs = vec![
            cfd_schematics::topology::presets::venturi_serpentine_series_spec(
                "series-option2",
                2.0e-3,
                0.7e-3,
                1.0e-3,
                1.8e-3,
                4,
                8.0e-3,
            ),
        ];
        let candidates =
            build_blueprint_candidates_from_specs(&specs, &operating_point).expect("candidates");

        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].treatment_channel_ids().len(), 5);
        assert_eq!(
            candidates[0].inferred_goal(),
            crate::OptimizationGoal::SelectiveVenturiCavitation
        );
    }

    #[test]
    fn canonical_parallel_specs_build_blueprint_candidates_without_bridge_logic() {
        let operating_point = operating_point(2.4e-6, 32_000.0, 0.12);
        let specs = vec![
            cfd_schematics::topology::presets::parallel_microchannel_array_spec(
                "parallel-option1",
                12,
                16.0e-3,
                150e-6,
                60e-6,
            ),
        ];
        let candidates =
            build_blueprint_candidates_from_specs(&specs, &operating_point).expect("candidates");

        assert_eq!(candidates.len(), 1);
        assert_eq!(candidates[0].treatment_channel_ids().len(), 12);
        let topology = candidates[0].topology_spec().expect("topology metadata");
        assert!(topology.has_parallel_paths());
        assert_eq!(
            candidates[0].inferred_goal(),
            crate::OptimizationGoal::SelectiveAcousticResidenceSeparation
        );
    }

    #[test]
    fn fresh_and_serialized_blueprint_metrics_match() {
        let candidate =
            canonical_option2_candidate("option2", operating_point(2.4e-6, 32_000.0, 0.12));
        let fresh = crate::evaluate_goal(
            &candidate,
            crate::OptimizationGoal::SelectiveVenturiCavitation,
        )
        .expect("fresh evaluation");
        let serialized = candidate.blueprint.to_json().expect("serialize");
        let restored_blueprint =
            cfd_schematics::NetworkBlueprint::from_json(&serialized).expect("roundtrip blueprint");
        let restored = crate::BlueprintCandidate::new(
            "option2-roundtrip",
            restored_blueprint,
            candidate.operating_point.clone(),
        );
        let roundtrip = crate::evaluate_goal(
            &restored,
            crate::OptimizationGoal::SelectiveVenturiCavitation,
        )
        .expect("roundtrip evaluation");

        assert!(
            (fresh.residence.treatment_residence_time_s
                - roundtrip.residence.treatment_residence_time_s)
                .abs()
                < 1.0e-12
        );
        assert!(
            (fresh.separation.separation_efficiency - roundtrip.separation.separation_efficiency)
                .abs()
                < 1.0e-12
        );
        assert!(
            (fresh.venturi.cavitation_selectivity_score
                - roundtrip.venturi.cavitation_selectivity_score)
                .abs()
                < 1.0e-12
        );
    }
}
