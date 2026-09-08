//! Derived topology queries computed from [`BlueprintTopologySpec`] structure.
//!
//! These methods replace the per-variant `match` dispatch formerly in
//! `cfd-optim::DesignTopology`. Every query is derived from the declarative
//! spec, not from an enum variant name, enabling the GA to compose arbitrary
//! topologies without extending an enum.

mod classification;
mod distribution;
mod naming;
mod venturi;

#[cfg(test)]
mod tests {
    use crate::topology::presets::enumerate_milestone12_topologies;
    use crate::TreatmentActuationMode;

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
}
