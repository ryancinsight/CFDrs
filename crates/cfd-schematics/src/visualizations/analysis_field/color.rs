//! Iris-to-schematic color boundary.

use iris::color::{ColorMap, NamedColorMap, Normalized};

use crate::error::{VisualizationError, VisualizationResult};
use crate::visualizations::traits::Color;

/// Map a scalar coordinate to the selected Iris color law.
///
/// Finite values outside `[0, 1]` are clamped. Iris evaluates the color law in
/// its native `f32` contract before the result is quantized to the schematic
/// backend's 8-bit channels.
///
/// # Errors
///
/// Returns [`VisualizationError::InvalidParameters`] when `value` is NaN or
/// infinite.
pub fn colorize(value: f64, map: NamedColorMap) -> VisualizationResult<Color> {
    if !value.is_finite() {
        return Err(VisualizationError::InvalidParameters {
            parameter: "normalized color coordinate".to_string(),
            value: value.to_string(),
            constraint: "value must be finite".to_string(),
        });
    }

    Ok(colorize_normalized(value.clamp(0.0, 1.0), map))
}

pub(super) fn colorize_normalized(value: f64, map: NamedColorMap) -> Color {
    let normalized = normalized_from_unit_interval(value);
    map.sample(normalized).into()
}

#[expect(
    clippy::cast_possible_truncation,
    reason = "Iris color laws intentionally evaluate normalized coordinates in f32"
)]
fn normalized_from_unit_interval(value: f64) -> Normalized {
    debug_assert!(value.is_finite() && (0.0..=1.0).contains(&value));
    Normalized::new(value as f32)
        .expect("invariant: finite clamped f64 remains normalized after conversion to f32")
}
