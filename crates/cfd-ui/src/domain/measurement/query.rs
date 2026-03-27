//! Measurement query types — state machine for multi-click measurement tools.

use crate::domain::scene::selection::SelectionTarget;

/// Describes what the user wants to measure.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum MeasurementQuery {
    /// Distance between two vertices/points (2 picks required).
    PointToPoint,
    /// Length of a selected edge (1 edge pick required).
    EdgeLength,
    /// Dihedral angle between two faces (2 face picks required).
    FaceToFaceAngle,
    /// Area of a selected face (1 face pick required).
    FaceArea,
    /// Volume of a closed mesh body (1 body pick required).
    MeshVolume,
}

impl MeasurementQuery {
    /// Number of picks required to complete this measurement.
    #[must_use]
    pub fn picks_required(self) -> usize {
        match self {
            Self::PointToPoint | Self::FaceToFaceAngle => 2,
            Self::EdgeLength | Self::FaceArea | Self::MeshVolume => 1,
        }
    }

    /// Unit string for the result of this measurement type.
    #[must_use]
    pub fn unit(self) -> &'static str {
        match self {
            Self::PointToPoint | Self::EdgeLength => "mm",
            Self::FaceToFaceAngle => "deg",
            Self::FaceArea => "mm^2",
            Self::MeshVolume => "mm^3",
        }
    }

    /// Display label for this measurement type.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::PointToPoint => "Point-to-Point Distance",
            Self::EdgeLength => "Edge Length",
            Self::FaceToFaceAngle => "Face Angle",
            Self::FaceArea => "Face Area",
            Self::MeshVolume => "Mesh Volume",
        }
    }
}

/// State machine for multi-click measurement tools.
#[derive(Clone, Debug, Default)]
pub enum MeasureToolState {
    /// No measurement in progress.
    #[default]
    Idle,
    /// Waiting for the first pick.
    AwaitingFirstPick(MeasurementQuery),
    /// First pick received, waiting for second.
    AwaitingSecondPick(MeasurementQuery, SelectionTarget),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_to_point_needs_two_picks() {
        assert_eq!(MeasurementQuery::PointToPoint.picks_required(), 2);
    }

    #[test]
    fn face_area_needs_one_pick() {
        assert_eq!(MeasurementQuery::FaceArea.picks_required(), 1);
    }

    #[test]
    fn units_are_correct() {
        assert_eq!(MeasurementQuery::PointToPoint.unit(), "mm");
        assert_eq!(MeasurementQuery::FaceToFaceAngle.unit(), "deg");
        assert_eq!(MeasurementQuery::FaceArea.unit(), "mm^2");
        assert_eq!(MeasurementQuery::MeshVolume.unit(), "mm^3");
    }
}
