//! View cube widget — interactive 3D orientation cube rendered in the viewport.
//!
//! Renders a small cube in the top-right corner of the viewport. Each face
//! is labeled with the corresponding view direction. Clicking a face/edge/corner
//! snaps the camera to the corresponding named engineering view.

use crate::domain::scene::named_views::NamedView;
use crate::domain::scene::view_cube::{ViewCubeRect, ViewCubeRegion, ViewCubeState};
use crate::infrastructure::gpu::overlay_buffer::{OverlayBuffer, OverlayTriangle};

/// Face colors for the view cube (slightly desaturated).
const FACE_COLORS: [[f32; 4]; 6] = [
    [0.35, 0.45, 0.55, 0.85], // Front
    [0.35, 0.45, 0.55, 0.85], // Back
    [0.45, 0.55, 0.35, 0.85], // Top
    [0.45, 0.55, 0.35, 0.85], // Bottom
    [0.55, 0.35, 0.45, 0.85], // Left
    [0.55, 0.35, 0.45, 0.85], // Right
];

/// Highlight tint applied when hovering a face.
const HOVER_TINT: [f32; 4] = [0.6, 0.7, 0.8, 0.95];

/// Generate view cube triangles for rendering.
///
/// Returns 12 triangles (2 per face) forming a unit cube centered at origin.
/// The `hovered` face gets a brighter color.
#[must_use]
pub fn view_cube_triangles(state: &ViewCubeState) -> Vec<OverlayTriangle> {
    let h = 0.5_f32;
    let faces: [([[f32; 3]; 4], usize); 6] = [
        // Front face (+Z)
        ([[- h, -h, h], [h, -h, h], [h, h, h], [-h, h, h]], 0),
        // Back face (-Z)
        ([[h, -h, -h], [-h, -h, -h], [-h, h, -h], [h, h, -h]], 1),
        // Top face (+Y)
        ([[-h, h, h], [h, h, h], [h, h, -h], [-h, h, -h]], 2),
        // Bottom face (-Y)
        ([[-h, -h, -h], [h, -h, -h], [h, -h, h], [-h, -h, h]], 3),
        // Left face (-X)
        ([[-h, -h, -h], [-h, -h, h], [-h, h, h], [-h, h, -h]], 4),
        // Right face (+X)
        ([[h, -h, h], [h, -h, -h], [h, h, -h], [h, h, h]], 5),
    ];

    let face_regions = [
        ViewCubeRegion::Front,
        ViewCubeRegion::Back,
        ViewCubeRegion::Top,
        ViewCubeRegion::Bottom,
        ViewCubeRegion::Left,
        ViewCubeRegion::Right,
    ];

    let mut tris = Vec::with_capacity(12);
    for (i, (verts, color_idx)) in faces.iter().enumerate() {
        let color = if state.hovered == Some(face_regions[i]) {
            HOVER_TINT
        } else {
            FACE_COLORS[*color_idx]
        };
        // Two triangles per quad face.
        tris.push(OverlayTriangle {
            v0: verts[0], v1: verts[1], v2: verts[2], color,
        });
        tris.push(OverlayTriangle {
            v0: verts[0], v1: verts[2], v2: verts[3], color,
        });
    }
    tris
}

/// Build GPU overlay buffer for the view cube.
pub fn build_view_cube_buffer(
    state: &ViewCubeState,
    device: &wgpu::Device,
) -> OverlayBuffer {
    OverlayBuffer::from_triangles(&view_cube_triangles(state), device)
}

/// Handle a click on the view cube. Returns the target `NamedView` if the
/// click landed on the cube, or `None` if it missed.
#[must_use]
pub fn handle_click(
    px: f32,
    py: f32,
    rect: &ViewCubeRect,
) -> Option<NamedView> {
    if !rect.contains(px, py) {
        return None;
    }
    let local_x = (px - rect.x) / rect.width;
    let local_y = (py - rect.y) / rect.height;
    let region = crate::domain::scene::view_cube::hit_test_face(local_x, local_y)?;
    Some(region.to_named_view())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_twelve_triangles() {
        let state = ViewCubeState::default();
        let tris = view_cube_triangles(&state);
        assert_eq!(tris.len(), 12);
    }

    #[test]
    fn hover_changes_color() {
        let mut state = ViewCubeState::default();
        let default_tris = view_cube_triangles(&state);
        state.hovered = Some(ViewCubeRegion::Front);
        let hover_tris = view_cube_triangles(&state);
        // First two triangles (front face) should have hover color.
        assert_ne!(default_tris[0].color, hover_tris[0].color);
    }
}
