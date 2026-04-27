//! Field legend panel — colorbar for simulation result visualization.

use crate::domain::colormap::{colorize, Colormap};

/// A colorbar legend entry (value -> color).
pub struct LegendEntry {
    /// Display label.
    pub label: String,
    /// RGBA color.
    pub color: [u8; 4],
}

/// Generate colorbar legend entries for a field range.
#[must_use]
pub fn generate_legend(
    min: f64,
    max: f64,
    steps: usize,
    colormap: Colormap,
    units: &str,
) -> Vec<LegendEntry> {
    (0..=steps)
        .map(|i| {
            let t = i as f64 / steps as f64;
            let value = min + t * (max - min);
            let color = colorize(t, colormap);
            LegendEntry {
                label: format!("{value:.3} {units}"),
                color: [color.r, color.g, color.b, color.a],
            }
        })
        .collect()
}
