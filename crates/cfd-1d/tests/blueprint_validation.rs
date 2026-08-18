#![expect(missing_docs, reason = "integration test crate exposes no public API")]

use aequitas::systems::si::quantities::Length;
use cfd_1d::validate_blueprint_for_1d_solve;
use cfd_schematics::geometry::generator::{
    create_primitive_selective_tree_geometry, PrimitiveSelectiveSplitKind,
    PrimitiveSelectiveTreeRequest,
};
fn selective_blueprint() -> cfd_schematics::NetworkBlueprint {
    create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
        name: "tri_tri".to_string(),
        box_dims_m: (Length::from_base(0.12776), Length::from_base(0.08547)),
        split_sequence: vec![
            PrimitiveSelectiveSplitKind::Tri,
            PrimitiveSelectiveSplitKind::Tri,
        ],
        main_width_m: Length::from_base(5.5e-3),
        throat_width_m: Length::from_base(45e-6),
        throat_length_m: Length::from_base(250e-6),
        channel_height_m: Length::from_base(1.0e-3),
        first_trifurcation_center_frac: 0.55,
        later_trifurcation_center_frac: 0.65,
        bifurcation_treatment_frac: 0.68,
        treatment_branch_venturi_enabled: true,
        treatment_branch_throat_count: 2,
        center_serpentine: None,
    })
}

#[test]
fn generated_selective_blueprint_passes_1d_graph_validation() {
    let blueprint = selective_blueprint();
    validate_blueprint_for_1d_solve(&blueprint)
        .expect("generated selective blueprint should validate");
}

#[test]
fn endpoint_inconsistent_channel_path_fails_1d_graph_validation() {
    let mut blueprint = selective_blueprint();
    blueprint.channels[0].path[0].0 += 2.0;

    let error = validate_blueprint_for_1d_solve(&blueprint)
        .expect_err("endpoint-inconsistent path must fail validation");
    assert!(
        error.to_string().contains("path start does not match"),
        "unexpected error: {error}"
    );
}
