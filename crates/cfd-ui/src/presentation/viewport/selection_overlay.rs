//! Selection overlay — renders highlight geometry for selected entities.
//!
//! Uses the `HighlightService` to generate overlay triangles (translucent
//! face highlights) and lines (edge highlights), then renders them via
//! the `OverlayPipeline` after the main mesh pass.

use crate::application::picking::highlight::HighlightService;
use crate::domain::scene::selection::{SelectionGranularity, SelectionSet};
use crate::infrastructure::gpu::overlay_buffer::OverlayBuffer;

/// Pre-built GPU buffers for the current selection highlight.
pub struct SelectionOverlay {
    /// Triangle overlay for face highlights (translucent blue).
    pub face_buffer: Option<OverlayBuffer>,
    /// Line overlay for edge highlights (bright white wireframe).
    pub edge_buffer: Option<OverlayBuffer>,
}

impl SelectionOverlay {
    /// Create an empty (no selection) overlay.
    #[must_use]
    pub fn empty() -> Self {
        Self {
            face_buffer: None,
            edge_buffer: None,
        }
    }

    /// Rebuild the overlay buffers from the current selection state.
    ///
    /// `face_positions` contains `[v0, v1, v2]` for each face of the selected node.
    pub fn rebuild(
        selection: &SelectionSet,
        granularity: SelectionGranularity,
        face_positions: &[[[f32; 3]; 3]],
        device: &wgpu::Device,
    ) -> Self {
        if selection.count() == 0 {
            return Self::empty();
        }

        let highlight_tris =
            HighlightService::highlight_faces(selection, granularity, face_positions);

        let face_buffer = if highlight_tris.is_empty() {
            None
        } else {
            Some(OverlayBuffer::from_triangles(&highlight_tris, device))
        };

        Self {
            face_buffer,
            edge_buffer: None, // Edge highlighting requires half-edge topology lookup.
        }
    }
}
