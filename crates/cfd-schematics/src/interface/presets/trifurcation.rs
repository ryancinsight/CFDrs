use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

const BLOOD_MU: f64 = 3.5e-3;

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
    bp.add_node(NodeSpec::new("junction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet_1", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("outlet_2", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("outlet_3", NodeKind::Outlet));

    bp.add_channel(ChannelSpec::new_pipe(
        "parent",
        "inlet",
        "junction",
        parent_length_m,
        parent_diameter_m,
        hp_resistance(parent_length_m, parent_diameter_m),
        0.0,
    ));

    for i in 1..=3 {
        bp.add_channel(ChannelSpec::new_pipe(
            format!("daughter_{i}"),
            "junction",
            format!("outlet_{i}"),
            daughter_length_m,
            daughter_diameter_m,
            hp_resistance(daughter_length_m, daughter_diameter_m),
            0.0,
        ));
    }

    bp
}

fn hp_resistance(length_m: f64, diameter_m: f64) -> f64 {
    128.0 * BLOOD_MU * length_m / (std::f64::consts::PI * diameter_m.powi(4))
}

/// Rectangular-channel symmetric trifurcation.
///
/// Topology: `inlet → parent → junction → (daughter_1, daughter_2, daughter_3) → outlets`.
/// All channels are rectangular with constant `height_m`.
#[must_use]
pub fn trifurcation_rect(
    name: impl Into<String>,
    parent_length_m: f64,
    daughter_length_m: f64,
    parent_width_m: f64,
    daughter_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("junction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet_1", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("outlet_2", NodeKind::Outlet));
    bp.add_node(NodeSpec::new("outlet_3", NodeKind::Outlet));

    bp.add_channel(ChannelSpec::new_pipe_rect(
        "parent",
        "inlet",
        "junction",
        parent_length_m,
        parent_width_m,
        height_m,
        shah_london_resistance(parent_width_m, height_m, parent_length_m, BLOOD_MU),
        0.0,
    ));

    for i in 1..=3 {
        bp.add_channel(ChannelSpec::new_pipe_rect(
            format!("daughter_{i}"),
            "junction",
            format!("outlet_{i}"),
            daughter_length_m,
            daughter_width_m,
            height_m,
            shah_london_resistance(daughter_width_m, height_m, daughter_length_m, BLOOD_MU),
            0.0,
        ));
    }

    bp
}

/// Shah (1978) Poiseuille-number resistance for a rectangular duct.
fn shah_london_resistance(width_m: f64, height_m: f64, length_m: f64, mu: f64) -> f64 {
    let alpha = height_m.min(width_m) / height_m.max(width_m);
    let po = 96.0
        * (1.0
            - 1.3553 * alpha
            + 1.9467 * alpha * alpha
            - 1.7012 * alpha.powi(3)
            + 0.9564 * alpha.powi(4)
            - 0.2537 * alpha.powi(5));
    let d_h = 2.0 * width_m * height_m / (width_m + height_m);
    let area = width_m * height_m;
    po * mu * length_m / (d_h * d_h * area)
}
