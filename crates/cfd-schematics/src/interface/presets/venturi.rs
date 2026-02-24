use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};

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

    bp.add_channel(
        ChannelSpec::new_pipe(
            "inlet_section",
            "inlet",
            "contraction",
            section_length,
            inlet_diameter_m,
            hp_resistance(section_length, inlet_diameter_m),
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
            hp_resistance(section_length, throat_diameter_m),
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
            hp_resistance(section_length, inlet_diameter_m),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

    bp
}

fn hp_resistance(length_m: f64, diameter_m: f64) -> f64 {
    128.0 * BLOOD_MU * length_m / (std::f64::consts::PI * diameter_m.powi(4))
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
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 5.0 * (inlet_width_m + throat_width_m) * 0.5;

    let mut bp = NetworkBlueprint::new(name);

    bp.add_node(NodeSpec::new("inlet", NodeKind::Inlet));
    bp.add_node(NodeSpec::new("contraction", NodeKind::Junction));
    bp.add_node(NodeSpec::new("throat", NodeKind::Junction));
    bp.add_node(NodeSpec::new("outlet", NodeKind::Outlet));

    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "inlet_section",
            "inlet",
            "contraction",
            l_taper,
            inlet_width_m,
            height_m,
            shah_london_resistance(inlet_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "throat_section",
            "contraction",
            "throat",
            l_throat,
            throat_width_m,
            height_m,
            shah_london_resistance(throat_width_m, height_m, l_throat, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::CancerTarget)),
    );
    bp.add_channel(
        ChannelSpec::new_pipe_rect(
            "diffuser_section",
            "throat",
            "outlet",
            l_taper,
            inlet_width_m,
            height_m,
            shah_london_resistance(inlet_width_m, height_m, l_taper, BLOOD_MU),
            0.0,
        )
        .with_metadata(TherapyZoneMetadata::new(TherapyZone::MixedFlow)),
    );

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
