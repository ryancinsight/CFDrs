//! Trifurcation-based composite presets (splits=3 wrappers around n_furcation).

use super::n_furcation::{n_furcation_serpentine_rect, n_furcation_venturi_rect};
use crate::domain::model::NetworkBlueprint;

/// Symmetric trifurcation with a venturi throat in each of 3 branches — closed loop.
#[must_use]
pub fn trifurcation_venturi_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> NetworkBlueprint {
    n_furcation_venturi_rect(
        name,
        3,
        trunk_length_m,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
    )
}

/// Symmetric trifurcation with a full serpentine in each of 3 arms — closed loop.
#[must_use]
pub fn trifurcation_serpentine_rect(
    name: impl Into<String>,
    trunk_length_m: f64,
    segments: usize,
    segment_length_m: f64,
    main_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    n_furcation_serpentine_rect(
        name,
        3,
        trunk_length_m,
        segments,
        segment_length_m,
        main_width_m,
        height_m,
    )
}
