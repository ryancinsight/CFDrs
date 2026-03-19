use crate::domain::{BlueprintCandidate, OperatingPoint};
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use cfd_schematics::topology::presets::{
    build_milestone12_blueprint, Milestone12StageBranchSpec, Milestone12StageLayout,
    Milestone12TopologyRequest,
};
use cfd_schematics::{
    BlueprintTopologySpec, BranchRole, SerpentineSpec, SplitKind, TreatmentActuationMode, VenturiPlacementMode,
};

pub(crate) fn canonical_option1_request() -> Milestone12TopologyRequest {
    let mut request = Milestone12TopologyRequest::new(
        "option1",
        "option1",
        vec![SplitKind::NFurcation(3), SplitKind::NFurcation(2)],
        6.0e-3,
        1.1e-3,
        22.0e-3,
        14.0e-3,
    );
    request.stage_layouts = vec![
        Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "wbc".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: 1.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 2.0e-3,
                },
                Milestone12StageBranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: 3.0e-3,
                },
            ],
        },
        Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(2),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 1.4e-3,
                },
                Milestone12StageBranchSpec {
                    label: "waste".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: 0.6e-3,
                },
            ],
        },
    ];
    request.center_serpentine = Some(SerpentineSpec {
        segments: 4,
        bend_radius_m: 1.3e-3,
        segment_length_m: 5.0e-3,
        wave_type: cfd_schematics::SerpentineWaveType::Sine,
    });
    request
}

pub(crate) fn canonical_option2_request() -> Milestone12TopologyRequest {
    let mut request = canonical_option1_request();
    request.topology_id = "option2".to_string();
    request.design_name = "option2".to_string();
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 2;
    request.venturi_throat_width_m = 0.4e-3;
    request.venturi_throat_length_m = 1.2e-3;
    request.venturi_placement_mode = VenturiPlacementMode::CurvaturePeakDeanNumber;
    request.venturi_target_channel_ids =
        vec![BlueprintTopologySpec::branch_channel_id("stage_1", "ctc")];
    request
}

pub(crate) fn canonical_option1_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
) -> BlueprintCandidate {
    let blueprint = build_milestone12_blueprint(&canonical_option1_request())
        .expect("canonical option1 fixture should build");
    BlueprintCandidate::new(id, blueprint, operating_point)
}

pub(crate) fn canonical_option2_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
) -> BlueprintCandidate {
    let blueprint = build_milestone12_blueprint(&canonical_option2_request())
        .expect("canonical option2 fixture should build");
    BlueprintCandidate::new(id, blueprint, operating_point)
}

pub(crate) fn stage0_venturi_candidate(
    id: impl Into<String>,
    operating_point: OperatingPoint,
    mode: VenturiPlacementMode,
) -> BlueprintCandidate {
    let mut request = Milestone12TopologyRequest::new(
        "venturi-fixture",
        "venturi-fixture",
        vec![SplitKind::NFurcation(3)],
        6.0e-3,
        1.1e-3,
        24.0e-3,
        14.0e-3,
    );
    request.stage_layouts = vec![Milestone12StageLayout {
        split_kind: SplitKind::NFurcation(3),
        branches: vec![
            Milestone12StageBranchSpec {
                label: "wbc".to_string(),
                role: BranchRole::Neutral,
                treatment_path: false,
                width_m: 1.0e-3,
            },
            Milestone12StageBranchSpec {
                label: "ctc".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                width_m: 2.0e-3,
            },
            Milestone12StageBranchSpec {
                label: "rbc".to_string(),
                role: BranchRole::RbcBypass,
                treatment_path: false,
                width_m: 3.0e-3,
            },
        ],
    }];
    request.center_serpentine = Some(SerpentineSpec {
        segments: 6,
        bend_radius_m: 1.1e-3,
        segment_length_m: 4.0e-3,
        wave_type: cfd_schematics::SerpentineWaveType::Sine,
    });
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 2;
    request.venturi_throat_width_m = 0.4e-3;
    request.venturi_throat_length_m = 1.2e-3;
    request.venturi_placement_mode = mode;
    request.venturi_target_channel_ids =
        vec![BlueprintTopologySpec::branch_channel_id("stage_0", "ctc")];

    let blueprint =
        build_milestone12_blueprint(&request).expect("stage0 venturi fixture should build");
    let treatment_lane_count = blueprint
        .channels
        .iter()
        .filter(|channel| channel.therapy_zone == Some(TherapyZone::CancerTarget))
        .count();
    assert!(
        treatment_lane_count >= 1,
        "fixture must preserve at least one treatment lane"
    );
    BlueprintCandidate::new(id, blueprint, operating_point)
}
