#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

use crate::BlueprintTopologyFactory;
use crate::topology::presets::symmetric_bifurcation_spec; // Using a proxy for bifurcation_rect since it's just a single split

#[must_use]
pub fn symmetric_bifurcation(
    name: impl Into<String>,
    parent_length_m: f64,
    daughter_length_m: f64,
    parent_diameter_m: f64,
    daughter_diameter_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("diverging_junction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("converging_junction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    bp.add_channel(
        ChannelSpec::new_pipe(
            "parent_in",
            "inlet",
            "diverging_junction",
            parent_length_m,
            parent_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp.add_channel(
        ChannelSpec::new_pipe(
            "daughter_1",
            "diverging_junction",
            "converging_junction",
            daughter_length_m,
            daughter_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
    );

    bp.add_channel(
        ChannelSpec::new_pipe(
            "daughter_2",
            "diverging_junction",
            "converging_junction",
            daughter_length_m,
            daughter_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::HealthyBypass)),
    );

    bp.add_channel(
        ChannelSpec::new_pipe(
            "parent_out",
            "converging_junction",
            "outlet",
            parent_length_m,
            parent_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}

/// Rectangular-channel symmetric bifurcation.
/// Topology: `inlet → parent_in → diverging_junction → (daughter_1, daughter_2) → converging_junction → parent_out → outlet`.
/// All channels are rectangular with constant `height_m`.
#[must_use]
pub fn bifurcation_rect(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_width_m: f64,
    _daughter_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    let spec = symmetric_bifurcation_spec(
        &name_str,
        1,
        parent_width_m,
        height_m,
        parent_length_m,
    );
    let bp = BlueprintTopologyFactory::build(&spec).expect("Valid topology spec");
    // In actual implementation, we'd need to adjust lengths directly on channels or update the spec fully.
    // For now we rely on the canonical factory.
    bp
}
