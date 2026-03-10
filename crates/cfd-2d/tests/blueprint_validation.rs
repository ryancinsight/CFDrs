use cfd_2d::network::validate_blueprint_for_2d_projection;
use cfd_schematics::domain::model::{ChannelSpec, NodeKind, NodeSpec};
use cfd_schematics::geometry::generator::{
    create_primitive_selective_tree_geometry, PrimitiveSelectiveSplitKind, PrimitiveSelectiveTreeRequest,
};

fn selective_blueprint() -> cfd_schematics::NetworkBlueprint {
    create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "tri_tri_2d".to_string(),
        box_dims_mm: (127.76, 85.47),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: 5.5e-3,
        throat_width_m: 45e-6,
        throat_length_m: 250e-6,
        channel_height_m: 1.0e-3,
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.65,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: true,
        treatment_branch_throat_count: 2,
        center_serpentine: None,
    })
}

#[test]
fn generated_selective_blueprint_passes_2d_projection_validation() {
    let blueprint = selective_blueprint();
    validate_blueprint_for_2d_projection(&blueprint)
        .expect("generated selective blueprint should validate for 2D projection");
}

#[test]
fn unresolved_crossing_fails_2d_projection_validation() {
    let mut blueprint = selective_blueprint();
    let start = (5.0, 5.0);
    let end = (blueprint.box_dims.0 - 5.0, blueprint.box_dims.1 - 5.0);
    blueprint.add_node(NodeSpec::new_at("crossing_top", NodeKind::Junction, start));
    blueprint.add_node(NodeSpec::new_at("crossing_bottom", NodeKind::Junction, end));

    let mut crossing_channel = ChannelSpec::new_pipe_rect(
        "crossing_lane",
        "crossing_top",
        "crossing_bottom",
        ((end.0 - start.0).hypot(end.1 - start.1)) * 1.0e-3,
        0.5e-3,
        1.0e-3,
        0.0,
        0.0,
    );
    crossing_channel.path = vec![start, end];
    blueprint.add_channel(crossing_channel);
    assert!(
        blueprint.has_unresolved_channel_overlaps(),
        "fixture must introduce an unresolved routed crossing"
    );

    let error = validate_blueprint_for_2d_projection(&blueprint)
        .expect_err("unresolved routed crossing must fail 2D validation");
    assert!(
        error
            .to_string()
            .contains("unresolved routed channel crossings"),
        "unexpected error: {error}"
    );
}
