//! Mixed split-family venturi bundles materialized through canonical topology specs.

use super::finalize_preset_blueprint;
use crate::domain::model::NetworkBlueprint;
use crate::domain::therapy_metadata::TherapyZone;
use crate::topology::presets::{parallel_path_spec, with_venturi};
use crate::topology::{
    ChannelRouteSpec, ParallelChannelSpec, ThroatGeometrySpec, TreatmentActuationMode,
    VenturiConfig, VenturiPlacementMode,
};
use crate::BlueprintTopologyFactory;

fn venturi_lane(
    channel_id: impl Into<String>,
    length_m: f64,
    width_m: f64,
    height_m: f64,
) -> ParallelChannelSpec {
    ParallelChannelSpec {
        channel_id: channel_id.into(),
        route: ChannelRouteSpec {
            length_m,
            width_m,
            height_m,
            serpentine: None,
            therapy_zone: TherapyZone::CancerTarget,
        },
    }
}

fn venturi_bundle_blueprint(
    name: String,
    inlet_width_m: f64,
    trunk_length_m: f64,
    channel_height_m: f64,
    throat_width_m: f64,
    throat_length_m: f64,
    channels: Vec<ParallelChannelSpec>,
) -> NetworkBlueprint {
    let spec = with_venturi(
        parallel_path_spec(
            &name,
            inlet_width_m,
            inlet_width_m,
            trunk_length_m,
            trunk_length_m,
            channels.clone(),
            TreatmentActuationMode::VenturiCavitation,
        ),
        VenturiConfig {
            target_channel_ids: channels
                .iter()
                .map(|channel| channel.channel_id.clone())
                .collect(),
            serial_throat_count: 1,
            throat_geometry: ThroatGeometrySpec {
                throat_width_m,
                throat_height_m: channel_height_m,
                throat_length_m,
                inlet_width_m: 0.0,
                outlet_width_m: 0.0,
                convergent_half_angle_deg: 7.0,
                divergent_half_angle_deg: 7.0,
            },
            placement_mode: VenturiPlacementMode::StraightSegment,
        },
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&spec)
            .expect("valid canonical mixed-tree venturi topology"),
    )
}

fn primary_lane_id(primary_label: &str, label: &str) -> String {
    if label == primary_label {
        "throat_section".to_string()
    } else {
        format!("throat_section_{label}")
    }
}

fn label_products(alphabet: &[&str], depth: usize) -> Vec<String> {
    if depth == 0 {
        return vec![String::new()];
    }
    let suffixes = label_products(alphabet, depth - 1);
    let mut labels = Vec::with_capacity(alphabet.len().pow(depth as u32));
    for prefix in alphabet {
        for suffix in &suffixes {
            labels.push(format!("{prefix}{suffix}"));
        }
    }
    labels
}

fn width_scaled_by_mixed_trifurcation(label: &str, main_width_m: f64, center_frac: f64) -> f64 {
    let center_frac = center_frac.clamp(0.20, 0.70);
    let periph_frac = (1.0 - center_frac) * 0.5;
    let level_1 = if label.starts_with('B') {
        center_frac
    } else {
        periph_frac
    };
    let level_2 = 0.5;
    let level_3 = 0.5;
    main_width_m
        * level_1
        * if label.len() >= 2 { level_2 } else { 1.0 }
        * if label.len() >= 3 { level_3 } else { 1.0 }
}

#[must_use]
pub fn double_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 2.5 * (main_width_m + throat_width_m);
    let leaf_length_m = 2.0 * branch1_length_m + 2.0 * l_taper + l_throat;
    let channels = label_products(&["A", "B", "C"], 2)
        .into_iter()
        .map(|label| {
            venturi_lane(
                primary_lane_id("AA", &label),
                leaf_length_m,
                main_width_m,
                height_m,
            )
        })
        .collect();
    venturi_bundle_blueprint(
        name.into(),
        main_width_m,
        trunk_length_m,
        height_m,
        throat_width_m,
        l_throat,
        channels,
    )
}

#[must_use]
pub fn bifurcation_trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 2.5 * (main_width_m + throat_width_m);
    let leaf_length_m = 2.0 * branch1_length_m + 2.0 * l_taper + l_throat;
    let labels = ["LA", "LB", "LC", "RA", "RB", "RC"];
    let channels = labels
        .into_iter()
        .map(|label| {
            venturi_lane(
                primary_lane_id("LA", label),
                leaf_length_m,
                main_width_m,
                height_m,
            )
        })
        .collect();
    venturi_bundle_blueprint(
        name.into(),
        main_width_m,
        trunk_length_m,
        height_m,
        throat_width_m,
        l_throat,
        channels,
    )
}

#[must_use]
pub fn trifurcation_bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 2.5 * (main_width_m + throat_width_m);
    let leaf_length_m = 2.0 * branch1_length_m + 2.0 * l_taper + l_throat;
    let labels = ["AL", "AR", "BL", "BR", "CL", "CR"];
    let channels = labels
        .into_iter()
        .map(|label| {
            venturi_lane(
                primary_lane_id("AL", label),
                leaf_length_m,
                main_width_m,
                height_m,
            )
        })
        .collect();
    venturi_bundle_blueprint(
        name.into(),
        main_width_m,
        trunk_length_m,
        height_m,
        throat_width_m,
        l_throat,
        channels,
    )
}

#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn trifurcation_bifurcation_bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    branch1_length_m: f64,
    branch2_length_m: f64,
    main_width_m: f64,
    center_frac: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    let l_throat = throat_length_m.max(2.0 * throat_width_m);
    let l_taper = 2.5 * (main_width_m + throat_width_m);
    let leaf_length_m = 2.0 * (branch1_length_m + branch2_length_m) + 2.0 * l_taper + l_throat;
    let channels = label_products(&["A", "B", "C"], 1)
        .into_iter()
        .flat_map(|root| {
            label_products(&["L", "R"], 2)
                .into_iter()
                .map(move |suffix| format!("{root}{suffix}"))
        })
        .map(|label| {
            venturi_lane(
                primary_lane_id("ALL", &label),
                leaf_length_m,
                width_scaled_by_mixed_trifurcation(&label, main_width_m, center_frac),
                height_m,
            )
        })
        .collect();
    venturi_bundle_blueprint(
        name.into(),
        main_width_m,
        trunk_length_m,
        height_m,
        throat_width_m,
        l_throat,
        channels,
    )
}
