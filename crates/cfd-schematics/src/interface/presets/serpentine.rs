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

        bp.add_channel(ChannelSpec::new_pipe(
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

/// Rectangular-channel serpentine network (photolithographic millifluidic chips).
///
/// All segments share the same cross-section `width_m × height_m`.
/// Resistance per segment uses the Shah (1978) Poiseuille-number correction for
/// rectangular ducts.  Segment channels are named `segment_1`, `segment_2`, …
#[must_use]
pub fn serpentine_rect(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    height_m: f64,
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

        bp.add_channel(ChannelSpec::new_pipe_rect(
            format!("segment_{}", i + 1),
            from,
            to,
            segment_length_m,
            width_m,
            height_m,
            shah_london_resistance(width_m, height_m, segment_length_m, BLOOD_MU),
            0.0,
        ));
    }

    bp
}

/// Shah (1978) Poiseuille-number resistance for a rectangular duct.
/// Returns `R` such that `ΔP = R · Q`.
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
