//! Parallel-path topology presets.

use crate::domain::therapy_metadata::TherapyZone;

use super::super::model::{BlueprintTopologySpec, ParallelChannelSpec, TreatmentActuationMode};
use super::helpers::{parallel_channel, PLATE_HEIGHT_MM, PLATE_WIDTH_MM};

/// Create a canonical parallel-channel topology spec.
#[must_use]
pub fn parallel_path_spec(
    name: &str,
    inlet_width_m: f64,
    outlet_width_m: f64,
    trunk_length_m: f64,
    outlet_tail_length_m: f64,
    parallel_channels: Vec<ParallelChannelSpec>,
    treatment_mode: TreatmentActuationMode,
) -> BlueprintTopologySpec {
    BlueprintTopologySpec {
        topology_id: format!("{name}_topology"),
        design_name: name.to_string(),
        box_dims_mm: (PLATE_WIDTH_MM, PLATE_HEIGHT_MM),
        inlet_width_m,
        outlet_width_m,
        trunk_length_m,
        outlet_tail_length_m,
        series_channels: Vec::new(),
        parallel_channels,
        split_stages: Vec::new(),
        venturi_placements: Vec::new(),
        treatment_mode,
    }
}

/// Canonical parallel microchannel leukapheresis array.
#[must_use]
pub fn parallel_microchannel_array_spec(
    name: &str,
    n_channels: usize,
    channel_length_m: f64,
    channel_width_m: f64,
    channel_height_m: f64,
) -> BlueprintTopologySpec {
    let channel_count = n_channels.max(1);
    let parallel_channels = (0..channel_count)
        .map(|index| {
            parallel_channel(
                format!("ch_{index}"),
                channel_length_m,
                channel_width_m,
                channel_height_m,
                TherapyZone::CancerTarget,
                None,
            )
        })
        .collect();
    parallel_path_spec(
        name,
        channel_width_m,
        channel_width_m,
        channel_length_m,
        channel_length_m,
        parallel_channels,
        TreatmentActuationMode::UltrasoundOnly,
    )
}
