//! Keyboard action definitions.

/// Named actions that can be bound to keyboard shortcuts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Action {
    NewProject,
    OpenProject,
    SaveProject,
    Undo,
    Redo,
    Delete,
    SelectAll,
    FitView,
    ToggleWireframe,
}

/// Default keyboard shortcut for an action (expressed as a string label).
#[must_use]
pub fn default_shortcut(action: Action) -> &'static str {
    match action {
        Action::NewProject => "Ctrl+N",
        Action::OpenProject => "Ctrl+O",
        Action::SaveProject => "Ctrl+S",
        Action::Undo => "Ctrl+Z",
        Action::Redo => "Ctrl+Shift+Z",
        Action::Delete => "Delete",
        Action::SelectAll => "Ctrl+A",
        Action::FitView => "F",
        Action::ToggleWireframe => "Z",
    }
}
