use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

const BLOOD_MU: f64 = 3.5e-3;

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

    bp.add_channel(ChannelSpec::new(
        "inlet_section",
        "inlet",
        "contraction",
        section_length,
        inlet_diameter_m,
        hp_resistance(section_length, inlet_diameter_m),
        0.0,
    ));

    bp.add_channel(ChannelSpec::new(
        "throat_section",
        "contraction",
        "throat",
        section_length,
        throat_diameter_m,
        hp_resistance(section_length, throat_diameter_m),
        0.0,
    ));

    bp.add_channel(ChannelSpec::new(
        "diffuser_section",
        "throat",
        "outlet",
        section_length,
        inlet_diameter_m,
        hp_resistance(section_length, inlet_diameter_m),
        0.0,
    ));

    bp
}

fn hp_resistance(length_m: f64, diameter_m: f64) -> f64 {
    128.0 * BLOOD_MU * length_m / (std::f64::consts::PI * diameter_m.powi(4))
}
