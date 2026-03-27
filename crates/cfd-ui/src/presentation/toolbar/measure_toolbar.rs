//! Measurement toolbar — tool buttons for each measurement type.

use crate::domain::measurement::MeasurementQuery;

/// Available measurement tool buttons.
#[must_use]
pub fn measurement_tools() -> Vec<MeasureToolButton> {
    vec![
        MeasureToolButton {
            query: MeasurementQuery::PointToPoint,
            label: "Distance",
            tooltip: "Measure point-to-point distance",
            shortcut: "D",
        },
        MeasureToolButton {
            query: MeasurementQuery::EdgeLength,
            label: "Edge",
            tooltip: "Measure edge length",
            shortcut: "L",
        },
        MeasureToolButton {
            query: MeasurementQuery::FaceToFaceAngle,
            label: "Angle",
            tooltip: "Measure angle between two faces",
            shortcut: "A",
        },
        MeasureToolButton {
            query: MeasurementQuery::FaceArea,
            label: "Area",
            tooltip: "Measure face area",
            shortcut: "R",
        },
        MeasureToolButton {
            query: MeasurementQuery::MeshVolume,
            label: "Volume",
            tooltip: "Compute mesh volume",
            shortcut: "V",
        },
    ]
}

/// Data for a single measurement tool button.
#[derive(Clone, Debug)]
pub struct MeasureToolButton {
    /// Which measurement this button activates.
    pub query: MeasurementQuery,
    /// Button label text.
    pub label: &'static str,
    /// Hover tooltip.
    pub tooltip: &'static str,
    /// Keyboard shortcut hint.
    pub shortcut: &'static str,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn five_measurement_tools() {
        assert_eq!(measurement_tools().len(), 5);
    }
}
