use super::{BlueprintTopologyFactory, BlueprintTopologyMutation};
use crate::domain::therapy_metadata::TherapyZone;
use crate::topology::presets::enumerate_milestone12_topologies;
use crate::{
    BlueprintTopologySpec, SerpentineSpec, SplitKind, TopologyOptimizationStage,
    TreatmentActuationMode, VenturiPlacementMode,
};
use aequitas::systems::si::quantities::{Angle, Length};

fn base_blueprint() -> crate::NetworkBlueprint {
    let request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "Tri-BASE")
        .expect("tri base request");
    crate::build_milestone12_blueprint(&request).expect("base blueprint")
}

fn treatment_venturi_channel_count(blueprint: &crate::NetworkBlueprint) -> usize {
    blueprint
        .channels
        .iter()
        .filter(|channel| {
            channel.venturi_geometry.is_some()
                && channel.therapy_zone == Some(TherapyZone::CancerTarget)
        })
        .count()
}

fn dean_placement(target_channel_id: String) -> crate::VenturiPlacementSpec {
    crate::VenturiPlacementSpec {
        placement_id: "test".to_string(),
        target_channel_id,
        serial_throat_count: 1,
        throat_geometry: crate::ThroatGeometrySpec {
            throat_width_m: Length::from_base(50.0e-6),
            throat_height_m: Length::from_base(1.0e-3),
            throat_length_m: Length::from_base(100.0e-6),
            inlet_width_m: Length::from_base(1.0e-3),
            outlet_width_m: Length::from_base(1.0e-3),
            convergent_half_angle: Angle::from_base(0.1),
            divergent_half_angle: Angle::from_base(0.1),
        },
        placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
    }
}

fn inserted_treatment_split_merge_blueprint(split_kind: SplitKind) -> crate::NetworkBlueprint {
    let blueprint = base_blueprint();
    let target_channel_id = blueprint
        .treatment_channel_ids()
        .into_iter()
        .next()
        .expect("treatment channel");
    let topology = blueprint.topology_spec().expect("topology");
    let route = topology
        .channel_route(&target_channel_id)
        .expect("treatment route");

    BlueprintTopologyFactory::mutate(
        &blueprint,
        BlueprintTopologyMutation::InsertTreatmentSplitMerge {
            target_channel_id,
            split_kind,
            treatment_serpentine: None,
            venturi_serial_throat_count: Some(1),
            venturi_throat_geometry: Some(crate::ThroatGeometrySpec {
                throat_width_m: Length::from_base(65.0e-6),
                throat_height_m: route.height_m,
                throat_length_m: Length::from_base(240.0e-6),
                inlet_width_m: route.width_m,
                outlet_width_m: route.width_m,
                convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
            }),
            venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .expect("split-merge mutation")
}

#[test]
fn treatment_channel_mutations_target_cancer_path_only() {
    let blueprint = base_blueprint();
    let target_channel_id = blueprint
        .treatment_channel_ids()
        .into_iter()
        .next()
        .expect("treatment channel");
    let topology = blueprint.topology_spec().expect("topology");
    let route = topology.channel_route(&target_channel_id).expect("route");

    let mutated = BlueprintTopologyFactory::mutate(
        &blueprint,
        BlueprintTopologyMutation::SetTreatmentChannelVenturi {
            target_channel_id: target_channel_id.clone(),
            serial_throat_count: 2,
            throat_geometry: crate::ThroatGeometrySpec {
                throat_width_m: Length::from_base(80.0e-6),
                throat_height_m: route.height_m,
                throat_length_m: Length::from_base(300.0e-6),
                inlet_width_m: route.width_m,
                outlet_width_m: route.width_m,
                convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
            },
            placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .expect("venturi mutation");

    assert!(
        mutated
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi)
    );
    assert!(mutated.channels.iter().any(|channel| channel.therapy_zone
        == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)));
}

#[test]
fn split_merge_insertion_preserves_geometry_authored_validation() {
    let mut request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "Quad-Y")
        .expect("quad mirror request");
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 1;
    let blueprint = crate::build_milestone12_blueprint(&request).expect("quad blueprint");
    let target_channel_id = blueprint
        .treatment_channel_ids()
        .into_iter()
        .next()
        .expect("treatment channel");

    let mutated = BlueprintTopologyFactory::mutate(
        &blueprint,
        BlueprintTopologyMutation::InsertTreatmentSplitMerge {
            target_channel_id,
            split_kind: SplitKind::NFurcation(3),
            treatment_serpentine: Some(SerpentineSpec {
                wave_type: crate::topology::SerpentineWaveType::Sine,
                segments: 4,
                bend_radius_m: Length::from_base(1.2e-3),
                segment_length_m: Length::from_base(4.0e-3),
            }),
            venturi_serial_throat_count: Some(2),
            venturi_throat_geometry: Some(crate::ThroatGeometrySpec {
                throat_width_m: Length::from_base(65.0e-6),
                throat_height_m: Length::from_base(1.0e-3),
                throat_length_m: Length::from_base(240.0e-6),
                inlet_width_m: Length::from_base(1.4e-3),
                outlet_width_m: Length::from_base(1.4e-3),
                convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
            }),
            venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .expect("split-merge mutation");

    assert!(mutated.is_geometry_authored());
    assert!(mutated.validate().is_ok());
    let topology = mutated.topology_spec().expect("mutated topology");
    assert!(topology.split_stages.len() >= 2);
    assert!(
        topology
            .venturi_placements
            .iter()
            .all(|placement| mutated.channels.iter().any(|channel| {
                (channel.id.as_str() == placement.target_channel_id
                    || channel
                        .id
                        .as_str()
                        .starts_with(&placement.target_channel_id))
                    && channel.venturi_geometry.is_some()
            })),
        "every declared venturi placement must materialize venturi geometry on a matching channel"
    );
}

#[test]
fn inserted_center_bifurcation_expands_venturis_across_both_child_channels() {
    let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(2));
    let topology = blueprint.topology_spec().expect("resolved topology");

    assert_eq!(treatment_venturi_channel_count(&blueprint), 2);
    assert_eq!(topology.venturi_placements.len(), 2);
    assert!(topology.venturi_placements.iter().all(|placement| {
        blueprint.channels.iter().any(|channel| {
            channel.id.as_str() == placement.target_channel_id
                && channel.therapy_zone == Some(TherapyZone::CancerTarget)
                && channel.venturi_geometry.is_some()
        })
    }));
}

#[test]
fn inserted_center_trifurcation_expands_venturis_across_all_child_channels() {
    let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(3));
    let topology = blueprint.topology_spec().expect("resolved topology");

    assert_eq!(treatment_venturi_channel_count(&blueprint), 3);
    assert_eq!(topology.venturi_placements.len(), 3);
    assert!(topology.venturi_placements.iter().all(|placement| {
        blueprint.channels.iter().any(|channel| {
            channel.id.as_str() == placement.target_channel_id
                && channel.therapy_zone == Some(TherapyZone::CancerTarget)
                && channel.venturi_geometry.is_some()
        })
    }));
}

#[test]
fn venturis_never_materialize_on_healthy_bypass_channels() {
    let blueprint = inserted_treatment_split_merge_blueprint(SplitKind::NFurcation(3));

    assert!(blueprint.channels.iter().all(|channel| {
        !(channel.therapy_zone == Some(TherapyZone::HealthyBypass)
            && channel.venturi_geometry.is_some())
    }));
}

#[test]
fn dean_site_prefers_target_route_serpentine_without_changing_values() {
    let mut blueprint = base_blueprint();
    let target_channel_id = blueprint
        .treatment_channel_ids()
        .into_iter()
        .next()
        .expect("treatment channel");
    blueprint
        .channels
        .iter_mut()
        .find(|channel| channel.therapy_zone == Some(TherapyZone::CancerTarget))
        .expect("materialized treatment channel")
        .id
        .0
        .clone_from(&target_channel_id);
    let mut topology = blueprint.topology_spec().expect("topology").clone();
    let expected_serpentine = SerpentineSpec {
        wave_type: crate::topology::SerpentineWaveType::Sine,
        segments: 3,
        bend_radius_m: Length::from_base(1.2e-3),
        segment_length_m: Length::from_base(4.0e-3),
    };
    for stage in &mut topology.split_stages {
        for branch in &mut stage.branches {
            if format!("{}_{}", stage.stage_id, branch.label) == target_channel_id {
                branch.route.serpentine = Some(expected_serpentine.clone());
            }
        }
    }
    let route = topology
        .channel_route(&target_channel_id)
        .expect("target route");
    let channel = blueprint
        .channels
        .iter()
        .find(|channel| channel.id.as_str() == target_channel_id)
        .expect("target channel");
    let area = channel.cross_section.area().into_base();
    let hydraulic_diameter = channel.cross_section.hydraulic_diameter().into_base();
    let average_velocity = 2.0e-9 / area;
    let reynolds = average_velocity * hydraulic_diameter / 4.0e-6;
    let serpentine = route.serpentine.as_ref().expect("target serpentine");
    let expected_radius = serpentine.bend_radius_m.into_base();
    let expected_arc_length = (serpentine.segments.max(1) as f64
        * (serpentine.segment_length_m.into_base().max(0.0)
            + std::f64::consts::PI * serpentine.bend_radius_m.into_base().max(0.0)))
    .max(channel.length_m.into_base());
    blueprint.topology = Some(topology);

    let estimate = BlueprintTopologyFactory::estimate_dean_site(
        &blueprint,
        &dean_placement(target_channel_id),
        2.0e-9,
        4.0e-6,
    )
    .expect("target route estimate");

    let expected_dean = reynolds * (hydraulic_diameter / (2.0 * expected_radius)).sqrt();
    assert_eq!(estimate.curvature_radius_m.into_base(), expected_radius);
    assert_eq!(estimate.arc_length_m.into_base(), expected_arc_length);
    assert_eq!(estimate.dean_number.into_base(), expected_dean);
}

#[test]
fn dean_site_uses_treatment_route_fallback_without_allocating_ids() {
    let mut blueprint = base_blueprint();
    let mut topology = blueprint.topology_spec().expect("topology").clone();
    let target_channel_id = topology
        .split_stages
        .iter()
        .flat_map(|stage| stage.branches.iter().map(move |branch| (stage, branch)))
        .find(|(_, branch)| !branch.treatment_path)
        .map(|(stage, branch)| format!("{}_{}", stage.stage_id, branch.label))
        .expect("healthy fallback route");
    let treatment_channel_id = topology
        .treatment_channel_ids()
        .into_iter()
        .next()
        .expect("treatment route");
    blueprint
        .channels
        .iter_mut()
        .find(|channel| channel.therapy_zone == Some(TherapyZone::HealthyBypass))
        .expect("materialized healthy fallback channel")
        .id
        .0
        .clone_from(&target_channel_id);

    let fallback_serpentine = SerpentineSpec {
        wave_type: crate::topology::SerpentineWaveType::Sine,
        segments: 2,
        bend_radius_m: Length::from_base(2.1e-3),
        segment_length_m: Length::from_base(3.0e-3),
    };
    for stage in &mut topology.split_stages {
        for branch in &mut stage.branches {
            if format!("{}_{}", stage.stage_id, branch.label) == treatment_channel_id {
                branch.route.serpentine = Some(fallback_serpentine.clone());
            }
            if format!("{}_{}", stage.stage_id, branch.label) == target_channel_id {
                branch.route.serpentine = None;
            }
        }
    }
    let route = topology
        .channel_route(&treatment_channel_id)
        .expect("fallback route");
    let channel = blueprint
        .channels
        .iter()
        .find(|channel| channel.id.as_str() == target_channel_id)
        .expect("target channel");
    let hydraulic_diameter = channel.cross_section.hydraulic_diameter().into_base();
    let area = channel.cross_section.area().into_base();
    let reynolds = (2.0e-9 / area) * hydraulic_diameter / 4.0e-6;
    let serpentine = route.serpentine.as_ref().expect("fallback serpentine");
    let expected_radius = serpentine.bend_radius_m.into_base();
    let expected_arc_length = (serpentine.segments.max(1) as f64
        * (serpentine.segment_length_m.into_base().max(0.0)
            + std::f64::consts::PI * serpentine.bend_radius_m.into_base().max(0.0)))
    .max(channel.length_m.into_base());
    blueprint.topology = Some(topology);

    let estimate = BlueprintTopologyFactory::estimate_dean_site(
        &blueprint,
        &dean_placement(target_channel_id),
        2.0e-9,
        4.0e-6,
    )
    .expect("treatment fallback estimate");

    let expected_dean = reynolds * (hydraulic_diameter / (2.0 * expected_radius)).sqrt();
    assert_eq!(estimate.curvature_radius_m.into_base(), expected_radius);
    assert_eq!(estimate.arc_length_m.into_base(), expected_arc_length);
    assert_eq!(estimate.dean_number.into_base(), expected_dean);
}
