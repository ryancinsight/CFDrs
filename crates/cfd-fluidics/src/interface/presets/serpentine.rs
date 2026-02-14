use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};

const BLOOD_MU: f64 = 3.5e-3;

#[must_use]
pub fn serpentine_chain(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    diameter_m: f64,
) -> NetworkBlueprint {
    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    for i in 1..segments {
        bp.add_node(NodeSpec::new(format!("turn_{i}"), NodeKind::Junction));
    }
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    for i in 0..segments {
        let from = if i == 0 {
            "inlet".to_string()
        } else {
            format!("turn_{i}")
        };

        let to = if i + 1 == segments {
            "outlet".to_string()
        } else {
            format!("turn_{}", i + 1)
        };

        bp.add_channel(ChannelSpec::new(
            format!("segment_{}", i + 1),
            from,
            to,
            segment_length_m,
            diameter_m,
            hp_resistance(segment_length_m, diameter_m),
            0.0,
        ));
    }

    bp
}

fn hp_resistance(length_m: f64, diameter_m: f64) -> f64 {
    128.0 * BLOOD_MU * length_m / (std::f64::consts::PI * diameter_m.powi(4))
}
