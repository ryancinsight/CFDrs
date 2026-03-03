//! Application workspace — the root layout containing all panels.

use crate::domain::document::project::ProjectDocument;
use crate::domain::document::history::CommandHistory;
use crate::domain::scene::selection::SelectionSet;
use crate::presentation::panels::console_panel::ConsoleState;
use crate::presentation::viewport::camera_controller::CameraController;
use crate::presentation::toolbar::mode_toolbar::EditMode;

/// Central application state holding all UI state and data.
pub struct WorkspaceState {
    /// The active project document.
    pub document: ProjectDocument,
    /// Command history for undo/redo.
    pub history: CommandHistory,
    /// Current object selection.
    pub selection: SelectionSet,
    /// Camera controller for the 3D viewport.
    pub camera_controller: CameraController,
    /// Console/log output.
    pub console: ConsoleState,
    /// Current editing mode.
    pub edit_mode: EditMode,
    /// Whether the application should exit.
    pub should_quit: bool,
}

impl WorkspaceState {
    /// Create a new workspace with default state.
    #[must_use]
    pub fn new() -> Self {
        Self {
            document: ProjectDocument::new(),
            history: CommandHistory::new(),
            selection: SelectionSet::new(),
            camera_controller: CameraController::new(),
            console: ConsoleState::new(1000),
            edit_mode: EditMode::default(),
            should_quit: false,
        }
    }
}

impl Default for WorkspaceState {
    fn default() -> Self {
        Self::new()
    }
}
