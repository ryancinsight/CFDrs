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
    let treatment_ids = seed.treatment_channel_ids();
    // Each treatment channel produces ~20 mutations (4 venturi + 1 serpentine
    // + 12 split-merge variants) plus 8 operating-point perturbations.
    let mut mutated = Vec::with_capacity(treatment_ids.len() * 20 + 8);

    for target_channel_id in treatment_ids {
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
        // Vary serpentine parameters: +1 segment, tighter bends, longer segments.
        // This produces genuinely different Dean numbers and residence times.
        let stronger_serpentine = Some(SerpentineSpec {
            segments: base_serpentine.segments + 2,
            bend_radius_m: base_serpentine.bend_radius_m * 0.8,
            segment_length_m: base_serpentine.segment_length_m * 1.2,
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

        // Venturi mutations with varied throat counts and widths.
        for &(vt_count, width_factor, label) in &[
            (2_u8, 1.0, "v2"),
            (3, 1.0, "v3"),
            (2, 0.7, "v2n"),  // narrower throat → lower σ
            (2, 1.3, "v2w"),  // wider throat → higher σ, less hemolysis
        ] {
            let mut varied_geom = venturi_geometry.clone();
            varied_geom.throat_width_m *= width_factor;
            let venturi_mutation = BlueprintTopologyFactory::mutate(
                &seed.blueprint,
                BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                    target_channel_id: target_channel_id.clone(),
                    serial_throat_count: vt_count,
                    throat_geometry: varied_geom,
                    placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
            )
            .map_err(OptimError::InvalidParameter)?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-{}-{}", seed.id, label, target_channel_id),
                venturi_mutation,
                seed.operating_point.clone(),
            ));
        }

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

    // Operating-point perturbation mutations: vary flow rate and pressure
    // to explore the physical design space around the seed topology.
    // These produce genuinely different scores because they change the
    // 1D CFD solve inputs (Re, ΔP, σ), unlike stacked serpentine
    // mutations which are nearly idempotent on the same channel.
    let base_q = seed.operating_point.flow_rate_m3_s;
    let base_p = seed.operating_point.inlet_gauge_pa;
    let base_ht = seed.operating_point.feed_hematocrit;
    for &q_factor in &[0.7, 0.85, 1.15, 1.3] {
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-qf{:.0}", seed.id, q_factor * 100.0),
            seed.blueprint.clone(),
            crate::domain::OperatingPoint {
                flow_rate_m3_s: base_q * q_factor,
                inlet_gauge_pa: base_p,
                feed_hematocrit: base_ht,
                patient_context: seed.operating_point.patient_context.clone(),
            },
        ));
    }
    for &p_factor in &[0.7, 0.85, 1.15, 1.3] {
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-pf{:.0}", seed.id, p_factor * 100.0),
            seed.blueprint.clone(),
            crate::domain::OperatingPoint {
                flow_rate_m3_s: base_q,
                inlet_gauge_pa: base_p * p_factor,
                feed_hematocrit: base_ht,
                patient_context: seed.operating_point.patient_context.clone(),
            },
        ));
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
