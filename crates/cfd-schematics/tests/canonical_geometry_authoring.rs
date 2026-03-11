use cfd_schematics::{
    build_milestone12_blueprint, build_milestone12_topology_spec, BlueprintTopologyFactory,
    BlueprintTopologySpec, BranchRole, Milestone12TopologyRequest, SplitKind,
    TreatmentActuationMode, VenturiPlacementMode,
};

fn selective_spec() -> BlueprintTopologySpec {
    let mut request = Milestone12TopologyRequest::new(
        "tri_tri_canonical",
        "tri_tri_canonical",
        vec![SplitKind::NFurcation(3), SplitKind::NFurcation(3)],
        5.5e-3,
        1.0e-3,
        20.0e-3,
        14.0e-3,
    );
    request.stage_layouts = vec![
        cfd_schematics::topology::presets::Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "wbc".to_string(),
                    role: BranchRole::WbcCollection,
                    treatment_path: false,
                    width_m: 1.0e-3,
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 2.0e-3,
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: 2.5e-3,
                },
            ],
        },
        cfd_schematics::topology::presets::Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "wbc".to_string(),
                    role: BranchRole::WbcCollection,
                    treatment_path: false,
                    width_m: 0.45e-3,
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: 1.3e-3,
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: 0.25e-3,
                },
            ],
        },
    ];
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 2;
    request.venturi_throat_width_m = 45e-6;
    request.venturi_throat_length_m = 250e-6;
    request.venturi_target_channel_ids =
        vec![BlueprintTopologySpec::branch_channel_id("stage_1", "ctc")];
    request.venturi_placement_mode = VenturiPlacementMode::StraightSegment;

    build_milestone12_topology_spec(&request)
}

#[test]
fn selective_factory_build_uses_canonical_geometry_authoring() {
    let blueprint = BlueprintTopologyFactory::build(&selective_spec())
        .expect("selective spec should build through canonical geometry path");

    assert!(
        blueprint.is_geometry_authored(),
        "selective factory builds must carry create_geometry provenance"
    );
    assert!(
        blueprint
            .channels
            .iter()
            .all(|channel| channel.path.len() >= 2),
        "report-grade channels must have explicit routed paths"
    );
    assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
    blueprint
        .validate()
        .expect("canonical selective blueprint should remain structurally valid");
}

#[test]
fn milestone12_blueprints_reject_missing_geometry_provenance() {
    let request = Milestone12TopologyRequest::new(
        "tri_manual_guard",
        "tri_manual_guard",
        vec![SplitKind::NFurcation(3)],
        6.0e-3,
        1.0e-3,
        10.0e-3,
        10.0e-3,
    );
    let mut blueprint =
        build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");
    blueprint.metadata = None;
    blueprint.geometry_authored = false;

    let error = blueprint
        .validate()
        .expect_err("manual provenance stripping must fail validation");
    assert!(error.contains("create_geometry"));
}
