use crate::domain::{BlueprintCandidate, OperatingPoint};
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::{
    BlueprintTopologyFactory, BlueprintTopologyMutation, BlueprintTopologySpec, BranchRole,
    BranchSpec, ChannelRouteSpec, SerpentineSpec, SplitKind, SplitStageSpec, ThroatGeometrySpec,
    TopologyOptimizationStage, TreatmentActuationMode, VenturiPlacementMode, VenturiPlacementSpec,
};

fn base_selective_spec(
    topology_id: &str,
    design_name: &str,
    split_stages: Vec<SplitStageSpec>,
    venturi_placements: Vec<VenturiPlacementSpec>,
    treatment_mode: TreatmentActuationMode,
) -> BlueprintTopologySpec {
    BlueprintTopologySpec {
        topology_id: topology_id.to_string(),
        design_name: design_name.to_string(),
        box_dims_mm: (127.76, 85.47),
        inlet_width_m: 6.0e-3,
        outlet_width_m: 4.0e-3,
        trunk_length_m: 22.0e-3,
        outlet_tail_length_m: 14.0e-3,
        series_channels: Vec::new(),
        parallel_channels: Vec::new(),
        split_stages,
        venturi_placements,
        treatment_mode,
    }
}

fn selective_trifurcation_stage(treatment_serpentine: Option<SerpentineSpec>) -> SplitStageSpec {
    SplitStageSpec {
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
                    height_m: 1.1e-3,
                    serpentine: None,
                    therapy_zone: TherapyZone::HealthyBypass,
                },
            },
            BranchSpec {
                label: "ctc".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                route: ChannelRouteSpec {
                    length_m: 28.0e-3,
                    width_m: 2.0e-3,
                    height_m: 1.1e-3,
                    serpentine: treatment_serpentine,
                    therapy_zone: TherapyZone::CancerTarget,
                },
            },
            BranchSpec {
                label: "rbc".to_string(),
                role: BranchRole::RbcBypass,
                treatment_path: false,
                route: ChannelRouteSpec {
                    length_m: 26.0e-3,
                    width_m: 3.0e-3,
                    height_m: 1.1e-3,
                    serpentine: None,
                    therapy_zone: TherapyZone::HealthyBypass,
                },
            },
        ],
    }
}

fn distal_treatment_stage() -> SplitStageSpec {
    SplitStageSpec {
        stage_id: "stage1".to_string(),
        split_kind: SplitKind::Bifurcation,
        branches: vec![
            BranchSpec {
                label: "ctc".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                route: ChannelRouteSpec {
                    length_m: 22.0e-3,
                    width_m: 1.4e-3,
                    height_m: 1.1e-3,
                    serpentine: Some(SerpentineSpec {
                        segments: 5,
                        bend_radius_m: 1.2e-3,
                        segment_length_m: 4.5e-3,
                    }),
                    therapy_zone: TherapyZone::CancerTarget,
                },
            },
            BranchSpec {
                label: "waste".to_string(),
                role: BranchRole::Neutral,
                treatment_path: false,
                route: ChannelRouteSpec {
                    length_m: 20.0e-3,
                    width_m: 0.6e-3,
                    height_m: 1.1e-3,
                    serpentine: None,
                    therapy_zone: TherapyZone::HealthyBypass,
                },
            },
        ],
    }
}

fn venturi_placement(target_channel_id: &str, mode: VenturiPlacementMode) -> VenturiPlacementSpec {
    VenturiPlacementSpec {
        placement_id: "p0".to_string(),
        target_channel_id: target_channel_id.to_string(),
        serial_throat_count: 2,
        throat_geometry: ThroatGeometrySpec {
            throat_width_m: 0.4e-3,
            throat_height_m: 1.1e-3,
            throat_length_m: 1.2e-3,
            inlet_width_m: 2.0e-3,
            outlet_width_m: 2.0e-3,
            convergent_half_angle_deg: 15.0,
            divergent_half_angle_deg: 9.0,
        },
        placement_mode: mode,
    }
}

pub(crate) fn canonical_option1_spec() -> BlueprintTopologySpec {
    base_selective_spec(
        "option1",
        "option1",
        vec![
            selective_trifurcation_stage(Some(SerpentineSpec {
                segments: 4,
                bend_radius_m: 1.3e-3,
                segment_length_m: 5.0e-3,
            })),
            distal_treatment_stage(),
        ],
        Vec::new(),
        TreatmentActuationMode::UltrasoundOnly,
    )
}

pub(crate) fn canonical_option1_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
) -> BlueprintCandidate {
    let spec = canonical_option1_spec();
    BlueprintCandidate::from_topology_spec(id, &spec, operating_point)
        .expect("canonical option1 fixture should build")
}

pub(crate) fn canonical_option2_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
) -> BlueprintCandidate {
    let option1_blueprint = BlueprintTopologyFactory::build(&canonical_option1_spec())
        .expect("canonical option1 fixture should build");
    let target_channel_id = BlueprintTopologySpec::branch_channel_id("stage1", "ctc");
    let blueprint = BlueprintTopologyFactory::mutate(
        &option1_blueprint,
        BlueprintTopologyMutation::UpdateVenturiConfiguration {
            placements: vec![venturi_placement(
                &target_channel_id,
                VenturiPlacementMode::CurvaturePeakDeanNumber,
            )],
            treatment_mode: TreatmentActuationMode::VenturiCavitation,
        },
        TopologyOptimizationStage::SelectiveVenturiCavitation,
    )
    .expect("canonical option2 fixture should mutate");
    BlueprintCandidate::new(id, blueprint, operating_point)
}

pub(crate) fn stage0_venturi_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
    mode: VenturiPlacementMode,
) -> BlueprintCandidate {
    let spec = base_selective_spec(
        "venturi-fixture",
        "venturi-fixture",
        vec![selective_trifurcation_stage(Some(SerpentineSpec {
            segments: 6,
            bend_radius_m: 1.1e-3,
            segment_length_m: 4.0e-3,
        }))],
        vec![venturi_placement(
            &BlueprintTopologySpec::branch_channel_id("stage0", "ctc"),
            mode,
        )],
        TreatmentActuationMode::VenturiCavitation,
    );
    BlueprintCandidate::from_topology_spec(id, &spec, operating_point)
        .expect("stage0 venturi fixture should build")
}
