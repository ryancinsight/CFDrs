//! Bifurcation-based composite presets (splits=2 wrappers around n_furcation).

use super::n_furcation::{n_furcation_serpentine_rect, n_furcation_venturi_rect};
use crate::domain::model::NetworkBlueprint;

/// Symmetric bifurcation with a venturi throat in each branch — closed loop.
#[must_use]
pub fn bifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    n_furcation_venturi_rect(
        name,
        2,
        trunk_length_m,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
    )
}

/// Symmetric bifurcation with a full serpentine in each arm — closed loop.
#[must_use]
pub fn bifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    n_furcation_serpentine_rect(
        name,
        2,
        trunk_length_m,
        segments,
        segment_length_m,
        main_width_m,
        height_m,
    )
}

/// Double-level symmetric bifurcation with 4 parallel serpentine arms — closed loop.
#[must_use]
pub fn double_bifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    // For double bifurcation we use n=4 splits (2×2 = 4 parallel arms)
    n_furcation_serpentine_rect(
        name,
        4,
        trunk_length_m,
        segments,
        segment_length_m,
        main_width_m,
        height_m,
    )
}
