//! Sketch toolbar — tool selector buttons for sketch mode.

use crate::presentation::viewport::sketch_interaction::SketchTool;

/// A button in the sketch toolbar.
#[derive(Clone, Debug)]
pub struct SketchToolButton {
    /// Which tool this button activates.
    pub tool: SketchTool,
    /// Button label.
    pub label: &'static str,
    /// Tooltip text.
    pub tooltip: &'static str,
    /// Keyboard shortcut hint.
    pub shortcut: &'static str,
}

/// Available sketch tool buttons.
#[must_use]
pub fn sketch_tools() -> Vec<SketchToolButton> {
    vec![
        SketchToolButton {
            tool: SketchTool::Select,
            label: "Select",
            tooltip: "Select/move sketch entities",
            shortcut: "S",
        },
        SketchToolButton {
            tool: SketchTool::Line,
            label: "Line",
            tooltip: "Draw a line segment",
            shortcut: "L",
        },
        SketchToolButton {
            tool: SketchTool::Arc3Point,
            label: "Arc",
            tooltip: "Draw a 3-point arc",
            shortcut: "A",
        },
        SketchToolButton {
            tool: SketchTool::Circle,
            label: "Circle",
            tooltip: "Draw a circle",
            shortcut: "C",
        },
        SketchToolButton {
            tool: SketchTool::Trim,
            label: "Trim",
            tooltip: "Trim to intersection",
            shortcut: "T",
        },
        SketchToolButton {
            tool: SketchTool::Offset,
            label: "Offset",
            tooltip: "Offset curve",
            shortcut: "O",
        },
        SketchToolButton {
            tool: SketchTool::Dimension,
            label: "Dimension",
            tooltip: "Add a driving dimension",
            shortcut: "D",
        },
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seven_sketch_tools() {
        assert_eq!(sketch_tools().len(), 7);
    }
}
