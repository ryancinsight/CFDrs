//! View cube — interactive orientation widget for the 3D viewport.
//!
//! The view cube is a 26-region interactive widget (6 faces + 12 edges +
//! 8 corners) rendered in a viewport corner. Each region maps to a named
//! engineering view that the camera snaps to on click.

use super::named_views::NamedView;

/// Identifies a clickable region on the view cube.
///
/// The 26 regions correspond to the faces, edges, and corners of a cube.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum ViewCubeRegion {
    // 6 faces
    Front,
    Back,
    Top,
    Bottom,
    Left,
    Right,
    // 12 edges (named by the two faces they connect)
    FrontTop,
    FrontBottom,
    FrontLeft,
    FrontRight,
    BackTop,
    BackBottom,
    BackLeft,
    BackRight,
    TopLeft,
    TopRight,
    BottomLeft,
    BottomRight,
    // 8 corners (named by three faces)
    FrontTopLeft,
    FrontTopRight,
    FrontBottomLeft,
    FrontBottomRight,
    BackTopLeft,
    BackTopRight,
    BackBottomLeft,
    BackBottomRight,
}

impl ViewCubeRegion {
    /// Map this region to the closest named engineering view.
    #[must_use]
    pub fn to_named_view(self) -> NamedView {
        match self {
            Self::Front => NamedView::Front,
            Self::Back => NamedView::Back,
            Self::Top => NamedView::Top,
            Self::Bottom => NamedView::Bottom,
            Self::Left => NamedView::Left,
            Self::Right => NamedView::Right,
            // Edge and corner regions snap to the nearest isometric view.
            Self::FrontTopRight | Self::FrontRight | Self::TopRight => {
                NamedView::IsoFrontTopRight
            }
            Self::FrontTopLeft | Self::FrontLeft | Self::TopLeft => {
                NamedView::IsoFrontTopLeft
            }
            Self::BackTopRight | Self::BackRight => NamedView::IsoBackTopRight,
            Self::BackTopLeft | Self::BackLeft => NamedView::IsoBackTopLeft,
            // Bottom edges/corners: snap to the face with the most visible area.
            Self::FrontTop => NamedView::Front,
            Self::FrontBottom | Self::FrontBottomLeft | Self::FrontBottomRight => {
                NamedView::Front
            }
            Self::BackTop | Self::BackBottom | Self::BackBottomLeft | Self::BackBottomRight => {
                NamedView::Back
            }
            Self::BottomLeft => NamedView::Left,
            Self::BottomRight => NamedView::Right,
        }
    }
}

/// Mutable state for the view cube widget.
#[derive(Clone, Debug, Default)]
pub struct ViewCubeState {
    /// The region currently hovered by the mouse, if any.
    pub hovered: Option<ViewCubeRegion>,
    /// Whether the camera is currently animating to a view-cube target.
    pub animating: bool,
}

/// Screen-space rectangle for the view cube (top-right corner of viewport).
#[derive(Clone, Copy, Debug)]
pub struct ViewCubeRect {
    /// Left edge in pixels from viewport left.
    pub x: f32,
    /// Top edge in pixels from viewport top.
    pub y: f32,
    /// Width in pixels.
    pub width: f32,
    /// Height in pixels.
    pub height: f32,
}

impl ViewCubeRect {
    /// Standard 120x120px rect in the top-right corner of a viewport.
    #[must_use]
    pub fn top_right(viewport_width: u32, viewport_height: u32) -> Self {
        let margin = 10.0;
        let size = 120.0;
        Self {
            x: viewport_width as f32 - size - margin,
            y: margin.min(viewport_height as f32 - size),
            width: size,
            height: size,
        }
    }

    /// Test whether a screen-space point is inside this rect.
    #[must_use]
    pub fn contains(&self, px: f32, py: f32) -> bool {
        px >= self.x && px <= self.x + self.width && py >= self.y && py <= self.y + self.height
    }
}

/// Hit-test a 2D point within the view cube rect to find the clicked region.
///
/// The cube is divided into a 3x3 grid on each visible face. The point
/// `(local_x, local_y)` is in `[0, 1]` normalized coordinates within the rect.
#[must_use]
pub fn hit_test_face(local_x: f32, local_y: f32) -> Option<ViewCubeRegion> {
    if !(0.0..=1.0).contains(&local_x) || !(0.0..=1.0).contains(&local_y) {
        return None;
    }

    let col = grid_cell(local_x);
    let row = grid_cell(local_y);

    // 3x3 grid: corners are (0,0), (0,2), (2,0), (2,2); edges are 1-cells.
    match (col, row) {
        (1, 1) => Some(ViewCubeRegion::Front), // center = face
        // Edges
        (1, 0) => Some(ViewCubeRegion::FrontTop),
        (1, 2) => Some(ViewCubeRegion::FrontBottom),
        (0, 1) => Some(ViewCubeRegion::FrontLeft),
        (2, 1) => Some(ViewCubeRegion::FrontRight),
        // Corners
        (0, 0) => Some(ViewCubeRegion::FrontTopLeft),
        (2, 0) => Some(ViewCubeRegion::FrontTopRight),
        (0, 2) => Some(ViewCubeRegion::FrontBottomLeft),
        (2, 2) => Some(ViewCubeRegion::FrontBottomRight),
        _ => None,
    }
}

/// Map a [0,1] coordinate to a grid cell index {0, 1, 2}.
fn grid_cell(t: f32) -> u8 {
    if t < EDGE_FRAC {
        0
    } else if t > 1.0 - EDGE_FRAC {
        2
    } else {
        1
    }
}

/// Fraction of cube face devoted to edge/corner hit regions.
const EDGE_FRAC: f32 = 0.25;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn center_hit_is_front_face() {
        assert_eq!(hit_test_face(0.5, 0.5), Some(ViewCubeRegion::Front));
    }

    #[test]
    fn corner_hit_is_front_top_left() {
        assert_eq!(hit_test_face(0.1, 0.1), Some(ViewCubeRegion::FrontTopLeft));
    }

    #[test]
    fn out_of_bounds_is_none() {
        assert_eq!(hit_test_face(-0.1, 0.5), None);
        assert_eq!(hit_test_face(0.5, 1.1), None);
    }

    #[test]
    fn top_right_rect_contains_center() {
        let rect = ViewCubeRect::top_right(800, 600);
        assert!(rect.contains(rect.x + 60.0, rect.y + 60.0));
        assert!(!rect.contains(0.0, 0.0));
    }

    #[test]
    fn face_regions_map_to_correct_views() {
        assert_eq!(ViewCubeRegion::Front.to_named_view(), NamedView::Front);
        assert_eq!(ViewCubeRegion::Back.to_named_view(), NamedView::Back);
        assert_eq!(ViewCubeRegion::Top.to_named_view(), NamedView::Top);
    }
}
