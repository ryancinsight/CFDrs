//! Render mode — shading and display mode types.

pub use crate::infrastructure::gpu::pipeline::ShadingMode;

/// Display mode for the viewport.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DisplayMode {
    /// Show the solid mesh with lighting.
    Shaded(ShadingMode),
    /// Show the mesh with field color overlay.
    FieldVisualization {
        /// Index into the simulation results field array.
        field_index: usize,
    },
}

impl Default for DisplayMode {
    fn default() -> Self {
        Self::Shaded(ShadingMode::Smooth)
    }
}
