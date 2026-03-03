//! Axis gizmo — coordinate axis indicator overlay.

/// Configuration for the axis gizmo displayed in the viewport corner.
#[derive(Clone, Debug)]
pub struct AxisGizmo {
    /// Size of the gizmo in pixels.
    pub size_px: u32,
    /// Whether the gizmo is visible.
    pub visible: bool,
    /// Position in the viewport (from top-left corner).
    pub offset_px: (u32, u32),
}

impl Default for AxisGizmo {
    fn default() -> Self {
        Self {
            size_px: 80,
            visible: true,
            offset_px: (20, 20),
        }
    }
}
