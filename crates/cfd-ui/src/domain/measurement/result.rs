//! Measurement result types — store completed measurement data.

/// Unique identifier for a measurement.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MeasurementId(pub u32);

/// The kind of measurement, with associated geometric data for rendering.
#[derive(Clone, Debug)]
pub enum MeasurementKind {
    /// Distance between two points.
    PointToPoint {
        start: [f64; 3],
        end: [f64; 3],
    },
    /// Length of an edge.
    EdgeLength {
        v0: [f64; 3],
        v1: [f64; 3],
    },
    /// Dihedral angle between two faces.
    FaceAngle {
        face_a_normal: [f64; 3],
        face_b_normal: [f64; 3],
        edge_midpoint: [f64; 3],
    },
    /// Area of a single face.
    FaceArea {
        face_center: [f64; 3],
    },
    /// Volume of a closed mesh.
    MeshVolume {
        mesh_center: [f64; 3],
    },
}

/// A completed measurement with its display data.
#[derive(Clone, Debug)]
pub struct MeasurementResult {
    /// Unique identifier.
    pub id: MeasurementId,
    /// What kind of measurement and its anchor geometry.
    pub kind: MeasurementKind,
    /// Scalar value of the measurement.
    pub value: f64,
    /// Unit string for display (e.g. "mm", "deg", "mm^2", "mm^3").
    pub unit: &'static str,
    /// Human-readable label.
    pub label: String,
}

impl MeasurementResult {
    /// Anchor points for overlay rendering (dimension line endpoints).
    #[must_use]
    pub fn anchor_points(&self) -> Vec<[f64; 3]> {
        match &self.kind {
            MeasurementKind::PointToPoint { start, end } => vec![*start, *end],
            MeasurementKind::EdgeLength { v0, v1 } => vec![*v0, *v1],
            MeasurementKind::FaceAngle { edge_midpoint, .. } => vec![*edge_midpoint],
            MeasurementKind::FaceArea { face_center } => vec![*face_center],
            MeasurementKind::MeshVolume { mesh_center } => vec![*mesh_center],
        }
    }

    /// Format the measurement value with unit for display.
    #[must_use]
    pub fn formatted_value(&self) -> String {
        format!("{:.4} {}", self.value, self.unit)
    }
}

/// Stores all active measurements for the session.
#[derive(Clone, Debug, Default)]
pub struct MeasurementStore {
    measurements: Vec<MeasurementResult>,
    next_id: u32,
}

impl MeasurementStore {
    /// Create an empty store.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Add a measurement and return its ID.
    pub fn add(
        &mut self,
        kind: MeasurementKind,
        value: f64,
        unit: &'static str,
        label: String,
    ) -> MeasurementId {
        let id = MeasurementId(self.next_id);
        self.next_id += 1;
        self.measurements.push(MeasurementResult {
            id,
            kind,
            value,
            unit,
            label,
        });
        id
    }

    /// Remove a measurement by ID.
    pub fn remove(&mut self, id: MeasurementId) {
        self.measurements.retain(|m| m.id != id);
    }

    /// Look up a measurement by ID.
    #[must_use]
    pub fn get(&self, id: MeasurementId) -> Option<&MeasurementResult> {
        self.measurements.iter().find(|m| m.id == id)
    }

    /// Iterate over all measurements.
    pub fn iter(&self) -> impl Iterator<Item = &MeasurementResult> {
        self.measurements.iter()
    }

    /// Number of measurements.
    #[must_use]
    pub fn count(&self) -> usize {
        self.measurements.len()
    }

    /// Clear all measurements.
    pub fn clear(&mut self) {
        self.measurements.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_and_retrieve_measurement() {
        let mut store = MeasurementStore::new();
        let id = store.add(
            MeasurementKind::PointToPoint {
                start: [0.0, 0.0, 0.0],
                end: [1.0, 0.0, 0.0],
            },
            1.0,
            "mm",
            "Distance".to_owned(),
        );
        assert_eq!(store.count(), 1);
        let m = store.get(id).expect("should exist");
        assert_eq!(m.value, 1.0);
        assert_eq!(m.unit, "mm");
    }

    #[test]
    fn remove_measurement() {
        let mut store = MeasurementStore::new();
        let id = store.add(
            MeasurementKind::FaceArea { face_center: [0.0; 3] },
            5.0,
            "mm^2",
            "Area".to_owned(),
        );
        store.remove(id);
        assert_eq!(store.count(), 0);
    }

    #[test]
    fn formatted_value() {
        let result = MeasurementResult {
            id: MeasurementId(0),
            kind: MeasurementKind::PointToPoint {
                start: [0.0; 3],
                end: [1.0, 0.0, 0.0],
            },
            value: 1.2345,
            unit: "mm",
            label: "Test".to_owned(),
        };
        assert_eq!(result.formatted_value(), "1.2345 mm");
    }
}
