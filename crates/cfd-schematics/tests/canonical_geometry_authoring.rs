//! Verifies canonical geometry authoring and its validation diagnostics.

use aequitas::systems::si::quantities::Length;
use cfd_schematics::geometry::generator::{
    create_selective_tree_geometry, SelectiveTreeRequest, SelectiveTreeTopology,
};
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
        Length::from_base(5.5e-3),
        Length::from_base(1.0e-3),
        Length::from_base(20.0e-3),
        Length::from_base(14.0e-3),
    );
    request.stage_layouts = vec![
        cfd_schematics::topology::presets::Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "wbc".to_string(),
                    role: BranchRole::WbcCollection,
                    treatment_path: false,
                    width_m: Length::from_base(1.0e-3),
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: Length::from_base(2.0e-3),
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: Length::from_base(2.5e-3),
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
                    width_m: Length::from_base(0.45e-3),
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "ctc".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: Length::from_base(1.3e-3),
                },
                cfd_schematics::topology::presets::Milestone12StageBranchSpec {
                    label: "rbc".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: Length::from_base(0.25e-3),
                },
            ],
        },
    ];
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 2;
    request.venturi_throat_width_m = Length::from_base(45e-6);
    request.venturi_throat_length_m = Length::from_base(250e-6);
    request.venturi_target_channel_ids =
        vec![BlueprintTopologySpec::branch_channel_id("stage_1", "ctc")];
    request.venturi_placement_mode = VenturiPlacementMode::StraightSegment;

    build_milestone12_topology_spec(&request)
}

#[test]
fn milestone12_request_quantities_propagate_to_typed_topology() {
    let request = Milestone12TopologyRequest::new(
        "typed-request",
        "typed-request",
        vec![SplitKind::NFurcation(2)],
        Length::from_base(5.5e-3),
        Length::from_base(1.0e-3),
        Length::from_base(20.0e-3),
        Length::from_base(14.0e-3),
    );

    let spec = build_milestone12_topology_spec(&request);
    assert_eq!(spec.box_dims_m.0, Length::from_base(127.76e-3));
    assert_eq!(spec.box_dims_m.1, Length::from_base(85.47e-3));
    assert_eq!(spec.inlet_width_m, request.inlet_width_m);
    assert_eq!(spec.trunk_length_m, request.branch_length_m);
    assert_eq!(spec.outlet_tail_length_m, request.outlet_tail_length_m);
    let (width_mm, height_mm) = request.box_dims_mm();
    let envelope_bound = 8.0 * f64::EPSILON * 128.0;
    assert!((width_mm - 127.76).abs() <= envelope_bound);
    assert!((height_mm - 85.47).abs() <= envelope_bound);
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
fn generic_selective_request_quantities_reach_geometry_builder() {
    let request = SelectiveTreeRequest {
        name: "typed-selective-request".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        trunk_length_m: Length::from_base(12.0e-3),
        branch_length_m: Length::from_base(10.0e-3),
        hybrid_branch_length_m: Length::from_base(8.0e-3),
        main_width_m: Length::from_base(1.2e-3),
        throat_width_m: Length::from_base(0.4e-3),
        throat_length_m: Length::from_base(3.0e-3),
        channel_height_m: Length::from_base(0.5e-3),
        topology: SelectiveTreeTopology::TriBiTriSelective {
            first_center_frac: 0.45,
            bi_treat_frac: 0.68,
            second_center_frac: 0.45,
        },
    };

    let blueprint = create_selective_tree_geometry(&request);

    let envelope_bound = 8.0 * f64::EPSILON * 128.0;
    assert!((blueprint.box_dims.0 - 127.76).abs() <= envelope_bound);
    assert!((blueprint.box_dims.1 - 85.47).abs() <= envelope_bound);
    assert!(blueprint.is_geometry_authored());
    assert!(blueprint
        .channels
        .iter()
        .all(|channel| channel.length_m.into_base().is_finite()));
}

#[test]
fn milestone12_blueprints_reject_missing_geometry_provenance() {
    let request = Milestone12TopologyRequest::new(
        "tri_manual_guard",
        "tri_manual_guard",
        vec![SplitKind::NFurcation(3)],
        Length::from_base(6.0e-3),
        Length::from_base(1.0e-3),
        Length::from_base(10.0e-3),
        Length::from_base(10.0e-3),
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
