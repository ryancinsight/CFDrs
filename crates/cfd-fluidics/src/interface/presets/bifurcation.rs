use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

const BLOOD_MU: f64 = 3.5e-3;

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
    bp.add_node(NodeSpec::new("junction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet_1", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("outlet_2", NodeKind::Outlet));

    bp.add_channel(ChannelSpec::new(
        "parent",
        "inlet",
        "junction",
        parent_length_m,
        parent_diameter_m,
        hp_resistance(parent_length_m, parent_diameter_m),
        0.0,
    ));

    bp.add_channel(ChannelSpec::new(
        "daughter_1",
        "junction",
        "outlet_1",
        daughter_length_m,
        daughter_diameter_m,
        hp_resistance(daughter_length_m, daughter_diameter_m),
        0.0,
    ));

    bp.add_channel(ChannelSpec::new(
        "daughter_2",
        "junction",
        "outlet_2",
        daughter_length_m,
        daughter_diameter_m,
        hp_resistance(daughter_length_m, daughter_diameter_m),
        0.0,
    ));

    bp
}

fn hp_resistance(length_m: f64, diameter_m: f64) -> f64 {
    128.0 * BLOOD_MU * length_m / (std::f64::consts::PI * diameter_m.powi(4))
}
