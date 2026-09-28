use crate::TreatmentActuationMode;
use crate::topology::presets::enumerate_milestone12_topologies;

#[test]
fn display_name_uses_actual_venturi_count_not_leaf_count() {
    let mut request = enumerate_milestone12_topologies()
        .into_iter()
        .find(|request| request.design_name == "PentaTriBi-BASE")
        .expect("PentaTriBi-BASE request");
    request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
    request.venturi_throat_count = 1;

    let blueprint = crate::build_milestone12_blueprint(&request).expect("venturi blueprint");
    let topology = blueprint.topology_spec().expect("resolved topology");

    let expected = format!(
        "{} + {}× Venturi",
        topology.stage_sequence_label(),
        topology.venturi_count()
    );
    let incorrect = format!(
        "{} + {}× Venturi",
        topology.stage_sequence_label(),
        topology.terminal_branch_count()
    );

    assert_eq!(topology.display_name(), expected);
    assert_ne!(topology.display_name(), incorrect);
}
