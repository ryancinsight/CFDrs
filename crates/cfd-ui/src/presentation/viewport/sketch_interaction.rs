//! Sketch interaction handler — mouse tool state machine for sketch mode.

use crate::domain::sketch::entity::EntityId;

/// Available sketch drawing tools.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum SketchTool {
    /// Select/move sketch entities.
    #[default]
    Select,
    /// Draw a line segment (two clicks).
    Line,
    /// Draw a 3-point arc.
    Arc3Point,
    /// Draw a circle (center + radius).
    Circle,
    /// Trim an entity at an intersection.
    Trim,
    /// Offset a curve by a distance.
    Offset,
    /// Place a driving dimension.
    Dimension,
}

/// State machine for sketch tool interactions.
#[derive(Clone, Debug, Default)]
pub enum SketchInteractionState {
    /// No active tool interaction.
    #[default]
    Idle,
    /// Line tool: first point placed, rubber-banding.
    LineFirstPoint(EntityId),
    /// Arc 3-point: first point placed.
    Arc3PointFirst(EntityId),
    /// Arc 3-point: two points placed.
    Arc3PointSecond(EntityId, EntityId),
    /// Circle: center placed, adjusting radius.
    CircleCenter(EntityId),
    /// Trim: hovering over candidate entity.
    TrimHover,
    /// Offset: entity selected, dragging distance.
    OffsetDragging(EntityId),
    /// Dimension: first reference picked.
    DimensionFirst(EntityId),
}

/// Handles mouse events for the active sketch tool.
pub struct SketchInteractionHandler {
    /// The currently active tool.
    pub tool: SketchTool,
    /// Current interaction state.
    pub state: SketchInteractionState,
    /// Snap radius in sketch-local units.
    pub snap_radius: f64,
}

impl Default for SketchInteractionHandler {
    fn default() -> Self {
        Self {
            tool: SketchTool::default(),
            state: SketchInteractionState::default(),
            snap_radius: 0.1,
        }
    }
}

impl SketchInteractionHandler {
    /// Create a new handler with default settings.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the active tool (resets state to Idle).
    pub fn set_tool(&mut self, tool: SketchTool) {
        self.tool = tool;
        self.state = SketchInteractionState::Idle;
    }

    /// Cancel the current tool interaction.
    pub fn cancel(&mut self) {
        self.state = SketchInteractionState::Idle;
    }

    /// Whether the handler is in the middle of a multi-click operation.
    #[must_use]
    pub fn is_active(&self) -> bool {
        !matches!(self.state, SketchInteractionState::Idle)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_tool_is_select() {
        let handler = SketchInteractionHandler::new();
        assert_eq!(handler.tool, SketchTool::Select);
        assert!(!handler.is_active());
    }

    #[test]
    fn set_tool_resets_state() {
        let mut handler = SketchInteractionHandler::new();
        handler.state = SketchInteractionState::LineFirstPoint(EntityId(0));
        handler.set_tool(SketchTool::Circle);
        assert_eq!(handler.tool, SketchTool::Circle);
        assert!(!handler.is_active());
    }
}
