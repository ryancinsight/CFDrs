//! Keyboard action definitions.

use crate::domain::measurement::{MeasurementId, MeasurementQuery};
use crate::domain::scene::named_views::NamedView;
use crate::domain::scene::selection::SelectionGranularity;

/// Named actions that can be bound to keyboard shortcuts.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum Action {
    // File operations
    NewProject,
    OpenProject,
    SaveProject,

    // Edit
    Undo,
    Redo,
    Delete,
    SelectAll,
    ClearSelection,

    // View
    FitView,
    ToggleWireframe,
    SnapToView(NamedView),
    ToggleViewCube,
    ToggleAxisIndicator,
    ToggleCameraMode,

    // Selection
    SetSelectionGranularity(SelectionGranularity),

    // Clipping
    AddClipPlane,
    RemoveClipPlane,

    // Measurement
    Measure(MeasurementQuery),
    ClearMeasurements,
    CopyMeasurement(MeasurementId),

    // Sketch
    EnterSketch,
    ExitSketch,
    SolveSketch,
}

/// Default keyboard shortcut for an action (expressed as a string label).
#[must_use]
pub fn default_shortcut(action: &Action) -> &'static str {
    match action {
        Action::NewProject => "Ctrl+N",
        Action::OpenProject => "Ctrl+O",
        Action::SaveProject => "Ctrl+S",
        Action::Undo => "Ctrl+Z",
        Action::Redo => "Ctrl+Shift+Z",
        Action::Delete => "Delete",
        Action::SelectAll => "Ctrl+A",
        Action::ClearSelection => "Escape",
        Action::FitView => "F",
        Action::ToggleWireframe => "Z",
        Action::SnapToView(NamedView::Front) => "1",
        Action::SnapToView(NamedView::Back) => "Ctrl+1",
        Action::SnapToView(NamedView::Top) => "7",
        Action::SnapToView(NamedView::Bottom) => "Ctrl+7",
        Action::SnapToView(NamedView::Left) => "3",
        Action::SnapToView(NamedView::Right) => "Ctrl+3",
        Action::SnapToView(NamedView::IsoFrontTopRight) => "0",
        Action::ToggleViewCube => "V",
        Action::ToggleAxisIndicator => "Shift+V",
        Action::ToggleCameraMode => "T",
        Action::SetSelectionGranularity(SelectionGranularity::Body) => "B",
        Action::SetSelectionGranularity(SelectionGranularity::Face) => "Shift+F",
        Action::SetSelectionGranularity(SelectionGranularity::Edge) => "E",
        Action::SetSelectionGranularity(SelectionGranularity::Vertex) => "Shift+E",
        _ => "",
    }
}
