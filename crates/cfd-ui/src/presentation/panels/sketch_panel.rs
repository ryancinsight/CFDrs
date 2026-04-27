//! Sketch panel — entity list, constraint list, and DOF indicator.

use crate::domain::sketch::dof::{DofAnalysis, DofStatus};
use crate::domain::sketch::sketch::Sketch;

/// Summary data for the sketch panel header.
#[derive(Clone, Debug)]
pub struct SketchPanelSummary {
    /// Number of entities.
    pub entity_count: usize,
    /// Number of constraints.
    pub constraint_count: usize,
    /// DOF analysis result.
    pub dof: DofAnalysis,
}

/// Build the panel summary from a sketch.
#[must_use]
pub fn sketch_panel_summary(sketch: &Sketch, dof: &DofAnalysis) -> SketchPanelSummary {
    SketchPanelSummary {
        entity_count: sketch.entity_count(),
        constraint_count: sketch.constraint_count(),
        dof: dof.clone(),
    }
}

/// Label for the DOF status indicator.
#[must_use]
pub fn dof_status_label(status: DofStatus) -> &'static str {
    match status {
        DofStatus::UnderConstrained => "Under-Constrained",
        DofStatus::FullyConstrained => "Fully Constrained",
        DofStatus::OverConstrained => "Over-Constrained",
    }
}

/// Color for the DOF indicator (CSS-style hex).
#[must_use]
pub fn dof_status_color(status: DofStatus) -> &'static str {
    match status {
        DofStatus::UnderConstrained => "#4488FF",
        DofStatus::FullyConstrained => "#44CC44",
        DofStatus::OverConstrained => "#FF4444",
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dof_labels_are_descriptive() {
        assert_eq!(
            dof_status_label(DofStatus::FullyConstrained),
            "Fully Constrained"
        );
        assert_eq!(
            dof_status_label(DofStatus::UnderConstrained),
            "Under-Constrained"
        );
        assert_eq!(
            dof_status_label(DofStatus::OverConstrained),
            "Over-Constrained"
        );
    }
}
