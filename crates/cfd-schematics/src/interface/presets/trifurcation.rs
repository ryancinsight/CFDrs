#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; nodes are created with NodeSpec::new_at().

use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

use crate::BlueprintTopologyFactory;
use crate::topology::presets::asymmetric_trifurcation_spec; // Using a proxy for trifurcation_rect since it's just a single split

#[must_use]
pub fn symmetric_trifurcation(
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

    bp.add_channel(ChannelSpec::new_pipe(
        "parent_in",
        "inlet",
        "diverging_junction",
        parent_length_m,
        parent_diameter_m,
        0.0, // resistance delegated to solver
        0.0,
    ));

    for i in 1..=3 {
        bp.add_channel(ChannelSpec::new_pipe(
            format!("daughter_{i}"),
            "diverging_junction",
            "converging_junction",
            daughter_length_m,
            daughter_diameter_m,
            0.0, // resistance delegated to solver
            0.0,
        ));
    }

    bp.add_channel(ChannelSpec::new_pipe(
        "parent_out",
        "converging_junction",
        "outlet",
        parent_length_m,
        parent_diameter_m,
        0.0, // resistance delegated to solver
        0.0,
    ));

    bp
}

/// Rectangular-channel symmetric trifurcation.
///
/// Topology: `inlet → parent → junction → (daughter_1, daughter_2, daughter_3) → outlets`.
/// All channels are rectangular with constant `height_m`.
#[must_use]
pub fn trifurcation_rect(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_width_m: f64,
    _daughter_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    // Assuming 1/3, 1/3, 1/3 for symmetric trifurcation.
    let center_frac = 1.0 / 3.0;
    let spec = asymmetric_trifurcation_spec(
        &name_str,
        1,
        parent_width_m,
        center_frac,
        height_m,
        parent_length_m,
    );
    BlueprintTopologyFactory::build(&spec).expect("Valid topology spec")
}
