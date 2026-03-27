//! Measurement panel — table of active measurements with values and actions.

use crate::domain::measurement::{MeasurementId, MeasurementResult, MeasurementStore};

/// Data for one row in the measurement results table.
#[derive(Clone, Debug)]
pub struct MeasureRow {
    /// Measurement identifier.
    pub id: MeasurementId,
    /// Measurement type label.
    pub kind_label: String,
    /// Formatted value with units.
    pub formatted_value: String,
    /// User label.
    pub label: String,
}

impl MeasureRow {
    /// Build a row from a measurement result.
    #[must_use]
    pub fn from_result(result: &MeasurementResult) -> Self {
        let kind_label = match &result.kind {
            crate::domain::measurement::MeasurementKind::PointToPoint { .. } => "Distance",
            crate::domain::measurement::MeasurementKind::EdgeLength { .. } => "Edge Length",
            crate::domain::measurement::MeasurementKind::FaceAngle { .. } => "Angle",
            crate::domain::measurement::MeasurementKind::FaceArea { .. } => "Area",
            crate::domain::measurement::MeasurementKind::MeshVolume { .. } => "Volume",
        };
        Self {
            id: result.id,
            kind_label: kind_label.to_owned(),
            formatted_value: result.formatted_value(),
            label: result.label.clone(),
        }
    }
}

/// Collect all measurement rows for the panel.
#[must_use]
pub fn measure_rows(store: &MeasurementStore) -> Vec<MeasureRow> {
    store.iter().map(MeasureRow::from_result).collect()
}

/// Format a measurement value for clipboard copying.
#[must_use]
pub fn format_for_clipboard(result: &MeasurementResult) -> String {
    format!("{}: {}", result.label, result.formatted_value())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::measurement::{MeasurementKind, MeasurementStore};

    #[test]
    fn measure_rows_from_store() {
        let mut store = MeasurementStore::new();
        store.add(
            MeasurementKind::PointToPoint {
                start: [0.0; 3],
                end: [1.0, 0.0, 0.0],
            },
            1.0,
            "mm",
            "Test".to_owned(),
        );
        let rows = measure_rows(&store);
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].kind_label, "Distance");
    }
}
