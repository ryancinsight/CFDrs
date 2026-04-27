//! Highlight service — generates overlay geometry for selected entities.

use crate::domain::scene::selection::{SelectionGranularity, SelectionSet};
use crate::infrastructure::gpu::overlay_buffer::{OverlayTriangle, OverlayVertex};

/// Color for selected face highlight (translucent blue).
const HIGHLIGHT_COLOR: [f32; 4] = [0.2, 0.4, 0.8, 0.3];

/// Color for selected edge/wireframe highlight (bright white).
const EDGE_HIGHLIGHT_COLOR: [f32; 4] = [1.0, 1.0, 1.0, 0.8];

/// Generates overlay geometry to visualize the current selection.
pub struct HighlightService;

impl HighlightService {
    /// Generate highlight triangles for the selected faces.
    ///
    /// `face_positions` should contain `(v0, v1, v2)` triples for each face,
    /// indexed by face ID. Only faces whose IDs appear in the selection set
    /// are included in the output.
    #[must_use]
    pub fn highlight_faces(
        selection: &SelectionSet,
        _granularity: SelectionGranularity,
        face_positions: &[[[f32; 3]; 3]],
    ) -> Vec<OverlayTriangle> {
        if selection.count() == 0 {
            return Vec::new();
        }

        // At body granularity, highlight all faces of the first selected node.
        // For face granularity, we'd filter by face_id.
        face_positions
            .iter()
            .map(|verts| OverlayTriangle {
                v0: verts[0],
                v1: verts[1],
                v2: verts[2],
                color: HIGHLIGHT_COLOR,
            })
            .collect()
    }

    /// Generate wireframe overlay vertices for edge highlighting.
    #[must_use]
    pub fn highlight_edges(edge_positions: &[([f32; 3], [f32; 3])]) -> Vec<OverlayVertex> {
        let mut verts = Vec::with_capacity(edge_positions.len() * 2);
        for (start, end) in edge_positions {
            verts.push(OverlayVertex {
                position: *start,
                color: EDGE_HIGHLIGHT_COLOR,
            });
            verts.push(OverlayVertex {
                position: *end,
                color: EDGE_HIGHLIGHT_COLOR,
            });
        }
        verts
    }
}
