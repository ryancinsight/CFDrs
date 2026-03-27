//! Clip plane panel — UI for managing section planes.
//!
//! Displays a list of 0-6 clip planes with enable/disable toggles,
//! normal/offset editors, and preset buttons for quick axis-aligned planes.

use crate::domain::clipping::{ClipPlane, ClipPlaneId, ClipPlaneSet, ClipPreset};

/// State for the clip plane panel UI.
#[derive(Clone, Debug, Default)]
pub struct ClipPanelState {
    /// Which plane is currently selected for editing, if any.
    pub selected: Option<ClipPlaneId>,
}

/// Data needed to render one row in the clip plane list.
#[derive(Clone, Debug)]
pub struct ClipPlaneRow {
    /// Plane identifier.
    pub id: ClipPlaneId,
    /// Display name.
    pub name: String,
    /// Whether clipping is active.
    pub enabled: bool,
    /// Whether the gizmo is visible.
    pub visible: bool,
    /// Normal vector components.
    pub normal: [f64; 3],
    /// Offset value.
    pub offset: f64,
}

impl ClipPlaneRow {
    /// Build a row from a clip plane reference.
    #[must_use]
    pub fn from_plane(plane: &ClipPlane) -> Self {
        Self {
            id: plane.id,
            name: plane.name.clone(),
            enabled: plane.enabled,
            visible: plane.visible,
            normal: [plane.normal.x, plane.normal.y, plane.normal.z],
            offset: plane.offset,
        }
    }
}

/// Collect all clip plane rows for the panel.
#[must_use]
pub fn clip_plane_rows(set: &ClipPlaneSet) -> Vec<ClipPlaneRow> {
    set.all_planes().map(ClipPlaneRow::from_plane).collect()
}

/// Available preset buttons.
#[must_use]
pub fn preset_buttons() -> Vec<(&'static str, ClipPreset)> {
    vec![
        ("XY", ClipPreset::Xy),
        ("XZ", ClipPreset::Xz),
        ("YZ", ClipPreset::Yz),
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_set_produces_no_rows() {
        let set = ClipPlaneSet::new();
        assert!(clip_plane_rows(&set).is_empty());
    }

    #[test]
    fn preset_buttons_count() {
        assert_eq!(preset_buttons().len(), 3);
    }
}
