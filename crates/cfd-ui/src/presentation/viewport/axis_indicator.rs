//! Axis indicator — RGB XYZ arrows rendered in a viewport corner (PyVista-style).
//!
//! Draws three perpendicular arrows (R=X, G=Y, B=Z) that rotate with the
//! camera orientation but ignore translation/zoom, giving a constant-size
//! orientation reference in the bottom-left corner.

use crate::infrastructure::gpu::overlay_buffer::{OverlayBuffer, OverlayLine};

/// Axis arrow colors: red=X, green=Y, blue=Z.
const X_COLOR: [f32; 4] = [0.9, 0.2, 0.2, 1.0];
const Y_COLOR: [f32; 4] = [0.2, 0.8, 0.2, 1.0];
const Z_COLOR: [f32; 4] = [0.2, 0.4, 0.9, 1.0];

/// Configuration for the axis indicator widget.
#[derive(Clone, Debug)]
pub struct AxisIndicatorConfig {
    /// Size of the indicator area in pixels.
    pub size_px: u32,
    /// Whether the indicator is visible.
    pub visible: bool,
    /// Offset from the bottom-left corner in pixels.
    pub offset_px: (u32, u32),
}

impl Default for AxisIndicatorConfig {
    fn default() -> Self {
        Self {
            size_px: 80,
            visible: true,
            offset_px: (20, 20),
        }
    }
}

/// Generate the axis indicator line geometry.
///
/// The arrows are unit-length lines from the origin along each axis.
/// The caller should transform these with the camera's rotation-only
/// matrix (no translation) before uploading to the GPU.
#[must_use]
pub fn axis_indicator_lines() -> Vec<OverlayLine> {
    let origin = [0.0, 0.0, 0.0];
    vec![
        OverlayLine {
            start: origin,
            end: [1.0, 0.0, 0.0],
            color: X_COLOR,
        },
        OverlayLine {
            start: origin,
            end: [0.0, 1.0, 0.0],
            color: Y_COLOR,
        },
        OverlayLine {
            start: origin,
            end: [0.0, 0.0, 1.0],
            color: Z_COLOR,
        },
    ]
}

/// Build GPU overlay buffer for the axis indicator.
#[must_use]
pub fn build_axis_buffer(device: &wgpu::Device) -> OverlayBuffer {
    OverlayBuffer::from_lines(&axis_indicator_lines(), device)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_three_axis_lines() {
        let lines = axis_indicator_lines();
        assert_eq!(lines.len(), 3);
    }

    #[test]
    fn default_config_is_visible() {
        let config = AxisIndicatorConfig::default();
        assert!(config.visible);
        assert_eq!(config.size_px, 80);
    }
}
