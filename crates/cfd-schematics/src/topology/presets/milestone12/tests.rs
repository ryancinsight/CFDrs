use crate::topology::model::{
    BranchRole, SplitKind, TreatmentActuationMode, VenturiPlacementMode,
};
use aequitas::systems::si::quantities::Length;

use super::super::sequence::MILESTONE12_SWEEP_SEQUENCES;
use super::catalog::MILESTONE12_MIRROR_VARIANTS;
use super::*;
use crate::topology::BlueprintTopologySpec;
use crate::topology::{
    BlueprintTopologyFactory, BlueprintTopologyMutation, TopologyOptimizationStage,
};

#[test]
fn default_stage_layouts_cover_bi_tri_quad_penta_roots() {
    let request = Milestone12PrimitiveSelectiveSpec::new(
        "m12",
        "m12",
        vec![
            SplitKind::NFurcation(2),
            SplitKind::NFurcation(3),
            SplitKind::NFurcation(4),
            SplitKind::NFurcation(5),
        ],
        Length::from_base(6.0e-3),
        Length::from_base(1.0e-3),
        Length::from_base(8.0e-3),
        Length::from_base(8.0e-3),
    );
    let layouts = milestone12_default_stage_layouts(&request);
    assert_eq!(layouts.len(), 4);
    assert_eq!(layouts[0].branches.len(), 2);
    assert_eq!(layouts[1].branches.len(), 3);
    assert_eq!(layouts[2].branches.len(), 4);
    assert_eq!(layouts[3].branches.len(), 5);
}

#[test]
fn explicit_even_way_treatment_lane_sets_are_preserved() {
    let mut request = Milestone12PrimitiveSelectiveSpec::new(
        "quad",
        "quad",
        vec![SplitKind::NFurcation(4)],
        Length::from_base(6.0e-3),
        Length::from_base(1.0e-3),
        Length::from_base(8.0e-3),
        Length::from_base(8.0e-3),
    );
    request.stage_layouts = vec![Milestone12StageLayout {
        split_kind: SplitKind::NFurcation(4),
        branches: vec![
            Milestone12StageBranchSpec {
                label: "arm_0".to_string(),
                role: BranchRole::Neutral,
                treatment_path: false,
                width_m: Length::from_base(1.0e-3),
            },
            Milestone12StageBranchSpec {
                label: "arm_1".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                width_m: Length::from_base(2.0e-3),
            },
            Milestone12StageBranchSpec {
                label: "arm_2".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                width_m: Length::from_base(2.0e-3),
            },
            Milestone12StageBranchSpec {
                label: "arm_3".to_string(),
                role: BranchRole::RbcBypass,
                treatment_path: false,
                width_m: Length::from_base(1.0e-3),
            },
        ],
    }];
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 1;
    request.venturi_target_channel_ids = vec![
        BlueprintTopologySpec::branch_channel_id("stage_0", "arm_1"),
        BlueprintTopologySpec::branch_channel_id("stage_0", "arm_2"),
    ];

    let spec = milestone12_primitive_selective_tree_spec(&request);
    assert_eq!(spec.treatment_channel_ids().len(), 2);
    assert_eq!(spec.venturi_placements.len(), 2);
}

#[test]
fn asymmetric_branch_widths_are_conserved() {
    let mut request = Milestone12PrimitiveSelectiveSpec::new(
        "tri",
        "tri",
        vec![SplitKind::NFurcation(3)],
        Length::from_base(6.0e-3),
        Length::from_base(1.0e-3),
        Length::from_base(8.0e-3),
        Length::from_base(8.0e-3),
    );
    request.stage_layouts = vec![Milestone12StageLayout {
        split_kind: SplitKind::NFurcation(3),
        branches: vec![
            Milestone12StageBranchSpec {
                label: "left".to_string(),
                role: BranchRole::Neutral,
                treatment_path: false,
                width_m: Length::from_base(1.0e-3),
            },
            Milestone12StageBranchSpec {
                label: "center".to_string(),
                role: BranchRole::Treatment,
                treatment_path: true,
                width_m: Length::from_base(2.0e-3),
            },
            Milestone12StageBranchSpec {
                label: "right".to_string(),
                role: BranchRole::RbcBypass,
                treatment_path: false,
                width_m: Length::from_base(3.0e-3),
            },
        ],
    }];

    let spec = milestone12_primitive_selective_tree_spec(&request);
    let total: f64 = spec.split_stages[0]
        .branches
        .iter()
        .map(|branch| branch.route.width_m.into_base())
        .sum();
    assert!((total - request.inlet_width_m.into_base()).abs() < 1.0e-12);
}

#[test]
fn enumerate_milestone12_topologies_covers_mirrored_root_catalog() {
    let catalog = enumerate_milestone12_topologies();
    assert_eq!(
        catalog.len(),
        MILESTONE12_SWEEP_SEQUENCES.len() * MILESTONE12_MIRROR_VARIANTS.len(),
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "Bi-BASE")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "Quad-Y")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "Penta-XY")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "QuadTriBi-BASE")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "PentaQuadBi-BASE")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "PentaQuadTri-XY")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "PentaTriBi-XY")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "QuadBi-X")
    );
    assert!(
        catalog
            .iter()
            .any(|request| request.design_name == "PentaTri-Y")
    );
}

#[test]
fn build_milestone12_blueprint_uses_canonical_geometry_authoring() {
    let request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "Tri-XY")
        .expect("mirrored Tri scaffold should exist");

    let blueprint =
        build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");

    assert!(blueprint.is_geometry_authored());
    assert!(
        blueprint
            .render_hints()
            .is_some_and(|hints| hints.mirror_x && hints.mirror_y)
    );
    blueprint
        .validate()
        .expect("Milestone 12 blueprint should validate");
}

#[test]
fn promote_option1_to_option2_rebuilds_geometry_authored_blueprint() {
    let mut request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "Quad-Y")
        .expect("mirrored Quad scaffold should exist");
    request.treatment_mode = TreatmentActuationMode::UltrasoundOnly;
    request.venturi_throat_count = 0;
    let option1 =
        build_milestone12_blueprint(&request).expect("Milestone 12 Option 1 should build");

    let promoted = promote_milestone12_option1_to_option2(
        &option1,
        1,
        VenturiPlacementMode::StraightSegment,
    )
    .expect("promotion should succeed");

    assert!(promoted.is_geometry_authored());
    assert!(
        promoted
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi)
    );
    assert!(
        promoted
            .render_hints()
            .is_some_and(|hints| !hints.mirror_x && hints.mirror_y)
    );
    assert_eq!(
        promoted
            .lineage()
            .and_then(|lineage| lineage.option1_source_blueprint.as_deref()),
        Some(option1.name.as_str())
    );
}

#[test]
fn selective_mutation_rejects_non_geometry_authored_blueprint() {
    let request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "Tri-BASE")
        .expect("Tri scaffold should exist");
    let mut blueprint =
        build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");
    blueprint.metadata = None;
    blueprint.geometry_authored = false;

    let error = BlueprintTopologyFactory::mutate(
        &blueprint,
        BlueprintTopologyMutation::UpdateBranchWidth {
            stage_id: "stage_0".to_string(),
            branch_label: "center".to_string(),
            new_width_m: 1.8e-3,
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .expect_err("non-canonical selective mutation must fail");
    assert!(
        error
            .to_string()
            .contains("create_geometry-authored provenance")
    );
}
