//! Phase 7 — End-to-end blueprint-native pipeline integration tests.
//!
//! These tests exercise the full path from `BlueprintTopologySpec` to
//! blueprint-native physics metrics, plus mesh generation from
//! `NetworkBlueprint`.

use cfd_optim::{evaluate_blueprint_candidate, BlueprintCandidate, OperatingPoint};
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::{
    BlueprintTopologyFactory, BlueprintTopologySpec, BranchRole, BranchSpec, ChannelRouteSpec,
    SerpentineSpec, SplitKind, SplitStageSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiPlacementMode, VenturiPlacementSpec,
};

fn selective_venturi_spec() -> BlueprintTopologySpec {
    BlueprintTopologySpec {
        topology_id: "pipeline-selective-venturi".to_string(),
        design_name: "pipeline-selective-venturi".to_string(),
        box_dims_mm: (127.76, 85.47),
        inlet_width_m: 4.0e-3,
        outlet_width_m: 3.0e-3,
        trunk_length_m: 20.0e-3,
        outlet_tail_length_m: 14.0e-3,
        split_stages: vec![SplitStageSpec {
            stage_id: "stage0".to_string(),
            split_kind: SplitKind::Trifurcation,
            branches: vec![
                BranchSpec {
                    label: "wbc".to_string(),
                    role: BranchRole::WbcCollection,
                    treatment_path: false,
                    route: ChannelRouteSpec {
                        length_m: 24.0e-3,
                        width_m: 1.0e-3,
                        height_m: 0.5e-3,
                        serpentine: None,
                        therapy_zone: TherapyZone::HealthyBypass,
                    },
                },
                BranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    route: ChannelRouteSpec {
                        length_m: 30.0e-3,
                        width_m: 2.0e-3,
                        height_m: 0.5e-3,
                        serpentine: Some(SerpentineSpec {
                            segments: 4,
                            bend_radius_m: 1.0e-3,
                            segment_length_m: 6.0e-3,
                        }),
                        therapy_zone: TherapyZone::CancerTarget,
                    },
                },
                BranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    route: ChannelRouteSpec {
                        length_m: 24.0e-3,
                        width_m: 1.0e-3,
                        height_m: 0.5e-3,
                        serpentine: None,
                        therapy_zone: TherapyZone::HealthyBypass,
                    },
                },
            ],
        }],
        venturi_placements: vec![VenturiPlacementSpec {
            placement_id: "vent-stage0".to_string(),
            target_channel_id: BlueprintTopologySpec::branch_channel_id("stage0", "ctc"),
            serial_throat_count: 1,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m: 50e-6,
                throat_height_m: 0.5e-3,
                throat_length_m: 100e-6,
                inlet_width_m: 2.0e-3,
                outlet_width_m: 2.0e-3,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 9.0,
            },
            placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
        }],
        treatment_mode: TreatmentActuationMode::VenturiCavitation,
    }
}

fn selective_venturi_candidate() -> BlueprintCandidate {
    BlueprintCandidate::new(
        "test-selective-venturi",
        BlueprintTopologyFactory::build(&selective_venturi_spec()).expect("blueprint must build"),
        OperatingPoint {
            flow_rate_m3_s: 1.0e-7,
            inlet_gauge_pa: 5_000.0,
            feed_hematocrit: 0.45,
            patient_context: None,
        },
    )
}

#[test]
fn selective_venturi_blueprint_candidate_computes_metrics() {
    let candidate = selective_venturi_candidate();
    let evaluation = evaluate_blueprint_candidate(&candidate)
        .expect("evaluate_blueprint_candidate must succeed");

    assert!(
        evaluation.venturi.placements[0]
            .cavitation_number
            .is_finite(),
        "cavitation_number must be finite, got {}",
        evaluation.venturi.placements[0].cavitation_number
    );
    assert!(
        evaluation.safety.pressure_drop_pa > 0.0,
        "pressure_drop_pa must be positive, got {}",
        evaluation.safety.pressure_drop_pa
    );
    assert!(
        evaluation.residence.treatment_residence_time_s > 0.0,
        "treatment_residence_time_s must be positive, got {}",
        evaluation.residence.treatment_residence_time_s
    );
}

#[test]
fn blueprint_json_round_trip_preserves_metrics() {
    let candidate = selective_venturi_candidate();
    let json = candidate
        .blueprint
        .to_json()
        .expect("blueprint JSON serialization must succeed");
    let roundtrip_blueprint =
        cfd_schematics::NetworkBlueprint::from_json(&json).expect("JSON roundtrip must succeed");

    assert_eq!(
        candidate.blueprint.nodes.len(),
        roundtrip_blueprint.nodes.len()
    );
    assert_eq!(
        candidate.blueprint.channels.len(),
        roundtrip_blueprint.channels.len()
    );
    assert_eq!(
        candidate.blueprint.treatment_channel_ids(),
        roundtrip_blueprint.treatment_channel_ids()
    );

    let fresh = evaluate_blueprint_candidate(&candidate).expect("fresh evaluation");
    let roundtrip = evaluate_blueprint_candidate(&BlueprintCandidate::new(
        "roundtrip-selective-venturi",
        roundtrip_blueprint,
        candidate.operating_point.clone(),
    ))
    .expect("roundtrip evaluation");

    assert!(
        (fresh.residence.treatment_residence_time_s
            - roundtrip.residence.treatment_residence_time_s)
            .abs()
            < 1.0e-12
    );
    assert!(
        (fresh.separation.separation_efficiency - roundtrip.separation.separation_efficiency).abs()
            < 1.0e-12
    );
}

#[test]
fn blueprint_mesh_pipeline_produces_watertight() {
    use cfd_mesh::application::pipeline::blueprint_mesh::{BlueprintMeshPipeline, PipelineConfig};
    use cfd_schematics::interface::presets::venturi_chain;

    let bp = venturi_chain("test", 0.030, 0.004, 0.002);
    let config = PipelineConfig::default();
    let mut output =
        BlueprintMeshPipeline::run(&bp, &config).expect("BlueprintMeshPipeline::run must succeed");

    assert!(
        output.fluid_mesh.is_watertight(),
        "Fluid mesh must be watertight"
    );
    assert!(
        output.fluid_mesh.face_count() > 0,
        "Fluid mesh must have faces, got {}",
        output.fluid_mesh.face_count()
    );
}
