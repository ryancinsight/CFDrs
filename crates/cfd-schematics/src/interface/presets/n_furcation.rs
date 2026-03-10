use super::finalize_preset_blueprint;
use crate::domain::model::NetworkBlueprint;
use crate::topology::presets;
use crate::topology::BlueprintTopologyFactory;
use crate::topology::BlueprintTopologySpec;

/// Generates a blueprint for a symmetric n-furcation.
#[must_use]
pub fn symmetric_n_furcation(
    name: &str,
    n_levels: usize,
    splits: usize,
    inlet_width_m: f64,
    channel_height_m: f64,
    trunk_length_m: f64,
) -> NetworkBlueprint {
    let spec = presets::symmetric_n_furcation_spec(
        name,
        n_levels,
        splits,
        inlet_width_m,
        channel_height_m,
        trunk_length_m,
    );

    n_furcation_rect(&spec)
}

/// Helper block constructing a general symmetric n-furcation network
/// from the given `BlueprintTopologySpec`.
#[must_use]
pub fn n_furcation_rect(spec: &BlueprintTopologySpec) -> NetworkBlueprint {
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(spec).expect("valid n_furcation topology spec"),
    )
}

// ── Backward-compatible bifurcation/trifurcation wrappers ──────────────────

/// Symmetric pipe-based bifurcation (1 inlet → 2 daughters → 1 outlet).
#[must_use]
pub fn symmetric_bifurcation(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_diameter_m: f64,
    _daughter_diameter_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    symmetric_n_furcation(
        &name_str,
        1,
        2,
        parent_diameter_m,
        parent_diameter_m,
        parent_length_m,
    )
}

/// Rectangular-channel symmetric bifurcation.
#[must_use]
pub fn bifurcation_rect(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_width_m: f64,
    _daughter_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    symmetric_n_furcation(&name_str, 1, 2, parent_width_m, height_m, parent_length_m)
}

/// Symmetric pipe-based trifurcation (1 inlet → 3 daughters → 1 outlet).
#[must_use]
pub fn symmetric_trifurcation(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_diameter_m: f64,
    _daughter_diameter_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    symmetric_n_furcation(
        &name_str,
        1,
        3,
        parent_diameter_m,
        parent_diameter_m,
        parent_length_m,
    )
}

/// Rectangular-channel symmetric trifurcation.
#[must_use]
pub fn trifurcation_rect(
    name: impl Into<String>,
    parent_length_m: f64,
    _daughter_length_m: f64,
    parent_width_m: f64,
    _daughter_width_m: f64,
    height_m: f64,
) -> NetworkBlueprint {
    let name_str = name.into();
    symmetric_n_furcation(&name_str, 1, 3, parent_width_m, height_m, parent_length_m)
}
