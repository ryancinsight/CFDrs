#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; all nodes are created with NodeSpec::new_at().

use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use crate::topology::presets::single_venturi_series_spec;
use crate::BlueprintTopologyFactory;



#[must_use]
pub fn venturi_chain(
    name: impl Into<String>,
    total_length_m: f64,
    inlet_diameter_m: f64,
    throat_diameter_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    let section_length = total_length_m / 3.0;

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    bp.add_channel(
        ChannelSpec::new_pipe(
            "inlet_section",
            "inlet",
            "contraction",
            section_length,
            inlet_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp.add_channel(
        ChannelSpec::new_pipe(
            "throat_section",
            "contraction",
            "throat",
            section_length,
            throat_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );

    bp.add_channel(
        ChannelSpec::new_pipe(
            "diffuser_section",
            "throat",
            "outlet",
            section_length,
            inlet_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}

/// Rectangular-channel planar venturi (photolithographic millifluidic chips).
///
/// The venturi narrows the channel **width** while keeping `height_m` constant —
/// the geometry used in PDMS/SU-8 chip fabrication.
///
/// Produces three channels:
/// - `inlet_section`:   `inlet_width_m × height_m`, converging approach
/// - `throat_section`:  `throat_width_m × height_m`, constant-area throat
/// - `diffuser_section`: `inlet_width_m × height_m`, expanding recovery
///
/// The converging/diverging section length is `5 × mean(inlet_width, throat_width)`.
/// Throat length is `max(throat_length_m, 2 × throat_width_m)`.
#[must_use]
pub fn venturi_rect(
    name: impl Into<String>,
    inlet_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let name = name.into();

    let spec = single_venturi_series_spec(
        &name,
        inlet_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
    );
    BlueprintTopologyFactory::build(&spec).expect("Valid topology spec")
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
            crate::TopologyOptimizationStage::SelectiveVenturiCavitation
        );
    }
}
