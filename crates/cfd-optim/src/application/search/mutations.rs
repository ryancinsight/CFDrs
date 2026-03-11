use crate::application::objectives::evaluate_goal;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use cfd_schematics::{
    promote_milestone12_option1_to_option2, BlueprintTopologyFactory, BlueprintTopologyMutation,
    SerpentineSpec, SplitKind, ThroatGeometrySpec, TopologyOptimizationStage,
    VenturiPlacementMode,
};

fn route_serpentine(route: &cfd_schematics::ChannelRouteSpec) -> SerpentineSpec {
    route.serpentine.clone().unwrap_or(SerpentineSpec {
        segments: 4,
        bend_radius_m: (route.width_m * 2.5).max(route.height_m),
        segment_length_m: (route.length_m / 4.0).max(route.width_m),
    })
}

fn route_throat_geometry(route: &cfd_schematics::ChannelRouteSpec) -> ThroatGeometrySpec {
    let inlet_width_m = route.width_m;
    let throat_width_m = (route.width_m * 0.22).clamp(35.0e-6, route.width_m * 0.85);
    ThroatGeometrySpec {
        throat_width_m,
        throat_height_m: route.height_m,
        throat_length_m: (route.length_m / 8.0).clamp(180.0e-6, route.length_m.max(220.0e-6)),
        inlet_width_m,
        outlet_width_m: inlet_width_m,
        convergent_half_angle_deg: 7.0,
        divergent_half_angle_deg: 7.0,
    }
}

pub fn seed_option2_candidates(candidates: &[BlueprintCandidate]) -> Vec<&BlueprintCandidate> {
    candidates
        .iter()
        .filter(|candidate| {
            candidate.inferred_goal()
                == OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity
                || candidate
                    .blueprint
                    .topology_spec()
                    .is_some_and(|spec| !spec.venturi_placements.is_empty())
        })
        .collect()
}

pub fn generate_ga_mutations(
    seed: &BlueprintCandidate,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let topology = seed.topology_spec()?;
    if topology.venturi_placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement requires canonical venturi-capable seeds, but candidate '{}' has no venturi placements",
            seed.id
        )));
    }
    let mut mutated = Vec::new();

    for target_channel_id in seed.treatment_channel_ids() {
        let route = topology
            .channel_route(&target_channel_id)
            .ok_or_else(|| {
                OptimError::InvalidParameter(format!(
                    "treatment channel '{}' not found in seed '{}'",
                    target_channel_id, seed.id
                ))
            })?
            .clone();
        let base_serpentine = route_serpentine(&route);
        let stronger_serpentine = Some(SerpentineSpec {
            segments: base_serpentine.segments + 1,
            bend_radius_m: base_serpentine.bend_radius_m,
            segment_length_m: base_serpentine.segment_length_m,
        });
        let venturi_geometry = route_throat_geometry(&route);

        let serpentine_mutation = BlueprintTopologyFactory::mutate(
            &seed.blueprint,
            BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
                target_channel_id: target_channel_id.clone(),
                serpentine: stronger_serpentine.clone(),
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .map_err(OptimError::InvalidParameter)?;
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-s-{}", seed.id, target_channel_id),
            serpentine_mutation,
            seed.operating_point.clone(),
        ));

        let venturi_mutation = BlueprintTopologyFactory::mutate(
            &seed.blueprint,
            BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                target_channel_id: target_channel_id.clone(),
                serial_throat_count: 2,
                throat_geometry: venturi_geometry.clone(),
                placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .map_err(OptimError::InvalidParameter)?;
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-v-{}", seed.id, target_channel_id),
            venturi_mutation,
            seed.operating_point.clone(),
        ));

        for split_kind in [
            SplitKind::NFurcation(2),
            SplitKind::NFurcation(3),
            SplitKind::NFurcation(4),
            SplitKind::NFurcation(5),
        ] {
            let split_merge = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: None,
                    venturi_serial_throat_count: None,
                    venturi_throat_geometry: None,
                    venturi_placement_mode: VenturiPlacementMode::StraightSegment,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-sm-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge,
                seed.operating_point.clone(),
            ));

            let split_merge_serpentine = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: stronger_serpentine.clone(),
                    venturi_serial_throat_count: None,
                    venturi_throat_geometry: None,
                    venturi_placement_mode: VenturiPlacementMode::StraightSegment,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-sms-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge_serpentine,
                seed.operating_point.clone(),
            ));

            let split_merge_venturi = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: stronger_serpentine.clone(),
                    venturi_serial_throat_count: Some(2),
                    venturi_throat_geometry: Some(venturi_geometry.clone()),
                    venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-smv-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge_venturi,
                seed.operating_point.clone(),
            ));
        }
    }

    Ok(mutated)
}

pub fn promote_option1_candidate_to_ga_seed(
    seed: &BlueprintCandidate,
) -> Result<BlueprintCandidate, OptimError> {
    let promoted = promote_milestone12_option1_to_option2(
        &seed.blueprint,
        1,
        VenturiPlacementMode::StraightSegment,
    )
    .map_err(OptimError::InvalidParameter)?;
    Ok(BlueprintCandidate::new(
        format!("{}-ga-promoted", seed.id),
        promoted,
        seed.operating_point.clone(),
    ))
}

pub fn build_milestone12_ga_seed_pair(
    candidates: &[BlueprintCandidate],
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let option1 = candidates
        .iter()
        .filter(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| !topology.has_venturi() && !topology.split_stages.is_empty())
        })
        .map(|candidate| {
            evaluate_goal(
                candidate,
                OptimizationGoal::AsymmetricSplitResidenceSeparation,
            )
            .map(|evaluation| (candidate, evaluation.score.unwrap_or(0.0)))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate.clone())
        .ok_or_else(|| {
            OptimError::InvalidParameter(
                "Milestone 12 GA seed selection found no Option 1 baseline candidate".to_string(),
            )
        })?;
    let option2 = candidates
        .iter()
        .filter(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| !topology.venturi_placements.is_empty())
        })
        .map(|candidate| {
            evaluate_goal(
                candidate,
                OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
            )
            .map(|evaluation| (candidate, evaluation.score.unwrap_or(0.0)))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate.clone())
        .ok_or_else(|| {
            OptimError::InvalidParameter(
                "Milestone 12 GA seed selection found no Option 2 baseline candidate".to_string(),
            )
        })?;

    Ok(vec![
        promote_option1_candidate_to_ga_seed(&option1)?,
        option2,
    ])
}
