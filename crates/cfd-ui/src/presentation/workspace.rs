//! Application workspace — the root layout containing all panels.

use crate::domain::clipping::ClipPlaneSet;
use crate::domain::document::history::CommandHistory;
use crate::domain::document::project::ProjectDocument;
use crate::domain::measurement::{MeasureToolState, MeasurementStore};
use crate::domain::scene::camera::animation::CameraAnimator;
use crate::domain::scene::selection::{SelectionGranularity, SelectionSet};
use crate::domain::scene::view_cube::ViewCubeState;
use crate::domain::sketch::SketchHandle;
use crate::presentation::panels::console_panel::ConsoleState;
use crate::presentation::toolbar::mode_toolbar::EditMode;
use crate::presentation::viewport::axis_indicator::AxisIndicatorConfig;
use crate::presentation::viewport::camera_controller::CameraController;
use crate::presentation::viewport::sketch_interaction::SketchInteractionHandler;

/// Central application state holding all UI state and data.
pub struct WorkspaceState {
    /// The active project document.
    pub document: ProjectDocument,
    /// Command history for undo/redo.
    pub history: CommandHistory,
    /// Current object selection.
    pub selection: SelectionSet,
    /// Selection granularity (body/face/edge/vertex).
    pub selection_granularity: SelectionGranularity,
    /// Camera controller for the 3D viewport.
    pub camera_controller: CameraController,
    /// Camera animation state (for view-cube snap transitions).
    pub camera_animator: CameraAnimator,
    /// View cube widget state.
    pub view_cube_state: ViewCubeState,
    /// Axis indicator configuration.
    pub axis_indicator: AxisIndicatorConfig,
    /// Active clip planes for section views.
    pub clip_planes: ClipPlaneSet,
    /// Active measurements.
    pub measurements: MeasurementStore,
    /// Measurement tool state machine.
    pub measure_tool_state: MeasureToolState,
    /// Active sketch handle (when in sketch edit mode).
    pub active_sketch: Option<SketchHandle>,
    /// Sketch tool interaction handler.
    pub sketch_interaction: SketchInteractionHandler,
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
            selection_granularity: SelectionGranularity::default(),
            camera_controller: CameraController::new(),
            camera_animator: CameraAnimator::new(),
            view_cube_state: ViewCubeState::default(),
            axis_indicator: AxisIndicatorConfig::default(),
            clip_planes: ClipPlaneSet::new(),
            measurements: MeasurementStore::new(),
            measure_tool_state: MeasureToolState::default(),
            active_sketch: None,
            sketch_interaction: SketchInteractionHandler::new(),
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
