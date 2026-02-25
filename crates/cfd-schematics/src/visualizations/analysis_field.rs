//! visualizations/analysis_field.rs — Typed CFD Analysis Field Overlay
//!
//! This module defines the domain-correct abstraction for overlaying CFD simulation
//! results onto 2D microfluidic schematics. Rather than accepting raw `Color` maps,
//! the renderer accepts a typed `AnalysisOverlay` that carries physical field data
//! and a colormap specification. Normalization and colorization are handled internally.
//!
//! # Invariants
//! - `AnalysisOverlay::edge_data` keys are channel IDs from `ChannelSystem::channels`.
//! - `AnalysisOverlay::node_data` keys are node IDs from `ChannelSystem::nodes`.
//! - `colorize` is a pure function: same inputs always produce the same `Color`.

use crate::visualizations::traits::Color;
use std::collections::HashMap;

/// Physical field quantity from a CFD simulation result
///
/// Each variant corresponds to a distinct physical observable that can be
/// mapped onto channels or nodes in a schematic visualization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AnalysisField {
    /// Static pressure [Pa] — node quantity, interpolated to edges for display
    Pressure,
    /// Wall shear stress [Pa] — edge quantity derived from velocity gradient at wall
    WallShearStress,
    /// Mean cross-sectional velocity [m/s] — edge quantity (Q / A)
    Velocity,
    /// Volumetric flow rate [m³/s] — edge quantity
    FlowRate,
    /// Apparent dynamic viscosity [Pa·s] — edge quantity (non-Newtonian fluids)
    Viscosity,
    /// User-defined field with a display label
    Custom(String),
}

impl AnalysisField {
    /// Human-readable label for axis/legend annotation
    #[must_use]
    pub fn label(&self) -> &str {
        match self {
            Self::Pressure => "Pressure [Pa]",
            Self::WallShearStress => "Wall Shear Stress [Pa]",
            Self::Velocity => "Velocity [m/s]",
            Self::FlowRate => "Flow Rate [m³/s]",
            Self::Viscosity => "Viscosity [Pa·s]",
            Self::Custom(s) => s.as_str(),
        }
    }
}

/// Colormap selection for field-to-color mapping
///
/// All colormaps map the normalized value `t ∈ [0, 1]` to an RGB color.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ColormapKind {
    /// Blue (low) → Red (high) — classic diverging heatmap
    #[default]
    BlueRed,
    /// Perceptually uniform sequential colormap (approximated)
    Viridis,
    /// Black (low) → White (high)
    Grayscale,
}

/// Map a normalized scalar `t ∈ [0, 1]` to a `Color` using the given colormap.
///
/// # Mathematical specification
/// - `BlueRed`: `R = round(255·t)`, `G = 0`, `B = round(255·(1−t))`
/// - `Viridis`: piecewise linear approximation of the Viridis LUT at 5 control points
/// - `Grayscale`: `R = G = B = round(255·t)`
///
/// Values outside `[0, 1]` are clamped.
#[must_use]
pub fn colorize(t: f64, colormap: ColormapKind) -> Color {
    let t = t.clamp(0.0, 1.0);
    match colormap {
        ColormapKind::BlueRed => Color::rgb(
            (255.0 * t) as u8,
            0,
            (255.0 * (1.0 - t)) as u8,
        ),
        ColormapKind::Viridis => {
            // 5-point piecewise linear approximation of Viridis
            // Control points: t=0.0, 0.25, 0.5, 0.75, 1.0
            const STOPS: &[(f64, u8, u8, u8)] = &[
                (0.00,  68,   1, 84),
                (0.25,  59,  82, 139),
                (0.50,  33, 145, 140),
                (0.75,  94, 201,  98),
                (1.00, 253, 231,  37),
            ];
            // Find surrounding stops
            let mut lo = STOPS[0];
            let mut hi = STOPS[STOPS.len() - 1];
            for window in STOPS.windows(2) {
                if t >= window[0].0 && t <= window[1].0 {
                    lo = window[0];
                    hi = window[1];
                    break;
                }
            }
            let span = hi.0 - lo.0;
            let s = if span > 0.0 { (t - lo.0) / span } else { 0.0 };
            let lerp = |a: u8, b: u8| -> u8 {
                (f64::from(a) + s * (f64::from(b) - f64::from(a))) as u8
            };
            Color::rgb(lerp(lo.1, hi.1), lerp(lo.2, hi.2), lerp(lo.3, hi.3))
        }
        ColormapKind::Grayscale => {
            let v = (255.0 * t) as u8;
            Color::rgb(v, v, v)
        }
    }
}

/// Complete overlay specification for rendering CFD results on a schematic
///
/// # Usage
/// ```rust,ignore
/// use cfd_schematics::visualizations::analysis_field::{AnalysisField, AnalysisOverlay, ColormapKind};
///
/// let overlay = AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::BlueRed)
///     .with_node_data(node_pressures)
///     .with_edge_data(edge_flow_rates);
/// renderer.render_analysis(&system, "output.png", &config, &overlay)?;
/// ```
#[derive(Debug, Clone)]
pub struct AnalysisOverlay {
    /// The physical field being visualized
    pub field: AnalysisField,
    /// Scalar values keyed by node ID (for node coloring, e.g. pressure at junctions)
    pub node_data: HashMap<usize, f64>,
    /// Scalar values keyed by channel ID (for edge coloring, e.g. flow rate, shear)
    pub edge_data: HashMap<usize, f64>,
    /// Colormap to use for normalization → color mapping
    pub colormap: ColormapKind,
}

impl AnalysisOverlay {
    /// Create a new overlay with the given field and colormap, no data yet.
    #[must_use]
    pub fn new(field: AnalysisField, colormap: ColormapKind) -> Self {
        Self {
            field,
            node_data: HashMap::new(),
            edge_data: HashMap::new(),
            colormap,
        }
    }

    /// Builder: attach node scalar data (e.g. pressure at each node ID).
    #[must_use]
    pub fn with_node_data(mut self, data: HashMap<usize, f64>) -> Self {
        self.node_data = data;
        self
    }

    /// Builder: attach edge scalar data (e.g. flow rate at each channel ID).
    #[must_use]
    pub fn with_edge_data(mut self, data: HashMap<usize, f64>) -> Self {
        self.edge_data = data;
        self
    }

    /// Compute the global min/max across all edge data values.
    ///
    /// Returns `(0.0, 1.0)` if no edge data is present (safe default).
    #[must_use]
    pub fn edge_range(&self) -> (f64, f64) {
        let values: Vec<f64> = self.edge_data.values().copied().collect();
        if values.is_empty() {
            return (0.0, 1.0);
        }
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        // Guard against degenerate case where all values are equal
        if (max - min).abs() < f64::EPSILON {
            (min - 1.0, max + 1.0)
        } else {
            (min, max)
        }
    }

    /// Compute the global min/max across all node data values.
    #[must_use]
    pub fn node_range(&self) -> (f64, f64) {
        let values: Vec<f64> = self.node_data.values().copied().collect();
        if values.is_empty() {
            return (0.0, 1.0);
        }
        let min = values.iter().copied().fold(f64::INFINITY, f64::min);
        let max = values.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        if (max - min).abs() < f64::EPSILON {
            (min - 1.0, max + 1.0)
        } else {
            (min, max)
        }
    }

    /// Get the `Color` for a given edge channel ID.
    ///
    /// Returns the default channel color (mid-gray) if the channel has no data.
    #[must_use]
    pub fn edge_color(&self, channel_id: usize) -> Option<Color> {
        let value = self.edge_data.get(&channel_id)?;
        let (min, max) = self.edge_range();
        let t = (value - min) / (max - min);
        Some(colorize(t, self.colormap))
    }

    /// Get the `Color` for a given node ID.
    #[must_use]
    pub fn node_color(&self, node_id: usize) -> Option<Color> {
        let value = self.node_data.get(&node_id)?;
        let (min, max) = self.node_range();
        let t = (value - min) / (max - min);
        Some(colorize(t, self.colormap))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colorize_blue_red_extremes() {
        let low = colorize(0.0, ColormapKind::BlueRed);
        assert_eq!(low, Color::rgb(0, 0, 255));

        let high = colorize(1.0, ColormapKind::BlueRed);
        assert_eq!(high, Color::rgb(255, 0, 0));
    }

    #[test]
    fn colorize_clamps_out_of_range() {
        let below = colorize(-1.0, ColormapKind::BlueRed);
        assert_eq!(below, colorize(0.0, ColormapKind::BlueRed));

        let above = colorize(2.0, ColormapKind::BlueRed);
        assert_eq!(above, colorize(1.0, ColormapKind::BlueRed));
    }

    #[test]
    fn colorize_grayscale_midpoint() {
        let mid = colorize(0.5, ColormapKind::Grayscale);
        assert_eq!(mid.r, mid.g);
        assert_eq!(mid.g, mid.b);
        assert_eq!(mid.r, 127);
    }

    #[test]
    fn overlay_edge_range_degenerate() {
        let mut overlay = AnalysisOverlay::new(AnalysisField::Pressure, ColormapKind::BlueRed);
        overlay.edge_data.insert(0, 100.0);
        overlay.edge_data.insert(1, 100.0);
        let (min, max) = overlay.edge_range();
        assert!(max > min, "degenerate range must be expanded");
    }

    #[test]
    fn overlay_edge_color_returns_none_for_missing_id() {
        let overlay = AnalysisOverlay::new(AnalysisField::FlowRate, ColormapKind::Viridis);
        assert!(overlay.edge_color(999).is_none());
    }
}
