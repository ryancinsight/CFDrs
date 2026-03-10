use super::finalize_preset_blueprint;
use crate::topology::presets::single_venturi_series_spec;
use crate::BlueprintTopologyFactory;

#[must_use]
pub fn venturi_chain(
    name: impl Into<String>,
    total_length_m: f64,
    inlet_diameter_m: f64,
    throat_diameter_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let throat_length_m = total_length_m / 3.0;
    let spec = single_venturi_series_spec(
        &name,
        inlet_diameter_m,
        throat_diameter_m,
        inlet_diameter_m,
        throat_length_m,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&spec).expect("valid venturi_chain topology spec"),
    )
}

/// Rectangular-channel planar venturi (photolithographic millifluidic chips).
#[must_use]
pub fn venturi_rect(
    name: impl Into<String>,
    inlet_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let spec = single_venturi_series_spec(
        &name,
        inlet_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&spec).expect("valid venturi_rect topology spec"),
    )
}

#[cfg(test)]
mod tests {
    use super::venturi_rect;

    #[test]
    fn venturi_rect_embeds_canonical_topology_metadata() {
        let blueprint = venturi_rect("rect-venturi", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3);
        let topology = blueprint
            .topology_spec()
            .expect("venturi_rect should attach topology metadata");
        let lineage = blueprint
            .lineage()
            .expect("venturi_rect should attach lineage metadata");

        assert_eq!(topology.design_name, "rect-venturi");
        assert_eq!(
            topology.treatment_channel_ids(),
            vec!["throat_section".to_string()]
        );
        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 1);
        assert_eq!(
            lineage.current_stage,
            crate::TopologyOptimizationStage::AsymmetricSplitVenturiCavitationSelectivity
        );
        assert!(blueprint.is_geometry_authored());
    }
}
