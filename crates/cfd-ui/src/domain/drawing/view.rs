//! Projected views — orthographic and section views for engineering drawings.

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

/// Standard and specialized view types for engineering drawings.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum ViewType {
    /// Standard orthographic: looking along -Z (front face).
    Front,
    /// Standard orthographic: looking along +Z (rear face).
    Back,
    /// Standard orthographic: looking along -Y (top-down).
    Top,
    /// Standard orthographic: looking along +Y (bottom-up).
    Bottom,
    /// Standard orthographic: looking along +X (left face).
    Left,
    /// Standard orthographic: looking along -X (right face).
    Right,
    /// Isometric projection (equal angles to all axes).
    Isometric,
    /// Section view cut by a plane.
    Section {
        /// Normal direction of the cutting plane.
        plane_normal: [f64; 3],
        /// Signed distance from origin along the normal.
        plane_offset: f64,
    },
    /// Detail view (magnified region).
    Detail {
        /// Center of the detail circle on the parent view (mm).
        center: [f64; 2],
        /// Radius of the detail circle on the parent view (mm).
        radius: f64,
        /// Magnification scale.
        detail_scale: f64,
    },
}

impl ViewType {
    /// The view direction vector (world space, looking from camera toward model).
    #[must_use]
    pub fn direction(&self) -> Vector3<f64> {
        match self {
            Self::Front => Vector3::new(0.0, 0.0, -1.0),
            Self::Back => Vector3::new(0.0, 0.0, 1.0),
            Self::Top => Vector3::new(0.0, -1.0, 0.0),
            Self::Bottom => Vector3::new(0.0, 1.0, 0.0),
            Self::Left => Vector3::new(1.0, 0.0, 0.0),
            Self::Right => Vector3::new(-1.0, 0.0, 0.0),
            Self::Isometric => {
                Vector3::new(-1.0, -1.0, -1.0).normalize()
            }
            Self::Section { plane_normal, .. } => {
                Vector3::new(plane_normal[0], plane_normal[1], plane_normal[2])
            }
            Self::Detail { .. } => Vector3::new(0.0, 0.0, -1.0),
        }
    }

    /// The up direction for this view.
    #[must_use]
    pub fn up(&self) -> Vector3<f64> {
        match self {
            Self::Top | Self::Bottom => Vector3::new(0.0, 0.0, -1.0),
            _ => Vector3::y(),
        }
    }
}

/// A projected view placed on a drawing sheet.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProjectedView {
    /// The type of projection.
    pub view_type: ViewType,
    /// Position of the view center on the sheet in millimeters from bottom-left.
    pub position_on_sheet_mm: [f64; 2],
    /// Drawing scale (e.g. 1.0 means 1:1, 0.5 means 1:2).
    pub scale: f64,
    /// Handle to the mesh being projected.
    pub mesh_index: usize,
    /// Optional label (e.g. "SECTION A-A", "DETAIL B").
    pub label: Option<String>,
}
