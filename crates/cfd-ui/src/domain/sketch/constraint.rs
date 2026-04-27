//! Sketch constraint types — geometric and dimensional constraints.

use super::entity::EntityId;

/// Unique identifier for a constraint within a sketch.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ConstraintId(pub u32);

/// Geometric and dimensional constraints between sketch entities.
#[derive(Clone, Debug)]
pub enum Constraint {
    /// A line is horizontal (start.y == end.y).
    Horizontal(EntityId),
    /// A line is vertical (start.x == end.x).
    Vertical(EntityId),
    /// Two points coincide (same position).
    Coincident(EntityId, EntityId),
    /// A line is tangent to an arc/circle at a shared point.
    Tangent(EntityId, EntityId),
    /// Two lines are perpendicular.
    Perpendicular(EntityId, EntityId),
    /// Two lines are parallel.
    Parallel(EntityId, EntityId),
    /// Two segments or radii have equal length.
    EqualLength(EntityId, EntityId),
    /// Two circles/arcs share a center point.
    Concentric(EntityId, EntityId),
    /// A point lies at the midpoint of a line.
    Midpoint { point: EntityId, line: EntityId },
    /// Two points are symmetric about a line.
    Symmetric {
        point_a: EntityId,
        point_b: EntityId,
        axis: EntityId,
    },
    /// A point is fixed in sketch space.
    Fixed(EntityId),
    /// Distance between two points (driving dimension).
    Distance {
        a: EntityId,
        b: EntityId,
        value: f64,
        driving: bool,
    },
    /// Angle between two lines.
    Angle {
        line_a: EntityId,
        line_b: EntityId,
        value_rad: f64,
        driving: bool,
    },
    /// Radius of a circle or arc.
    Radius {
        entity: EntityId,
        value: f64,
        driving: bool,
    },
}

impl Constraint {
    /// Number of scalar residual equations this constraint contributes.
    #[must_use]
    pub fn residual_count(&self) -> usize {
        match self {
            Self::Horizontal(_)
            | Self::Vertical(_)
            | Self::Tangent(_, _)
            | Self::Perpendicular(_, _)
            | Self::Parallel(_, _)
            | Self::EqualLength(_, _)
            | Self::Distance { .. }
            | Self::Angle { .. }
            | Self::Radius { .. } => 1,

            Self::Coincident(_, _)
            | Self::Concentric(_, _)
            | Self::Fixed(_)
            | Self::Midpoint { .. }
            | Self::Symmetric { .. } => 2,
        }
    }

    /// Entity IDs referenced by this constraint.
    #[must_use]
    pub fn referenced_entities(&self) -> Vec<EntityId> {
        match self {
            Self::Horizontal(e)
            | Self::Vertical(e)
            | Self::Fixed(e)
            | Self::Radius { entity: e, .. } => {
                vec![*e]
            }
            Self::Coincident(a, b)
            | Self::Tangent(a, b)
            | Self::Perpendicular(a, b)
            | Self::Parallel(a, b)
            | Self::EqualLength(a, b)
            | Self::Concentric(a, b)
            | Self::Distance { a, b, .. }
            | Self::Angle {
                line_a: a,
                line_b: b,
                ..
            } => {
                vec![*a, *b]
            }
            Self::Midpoint { point, line } => vec![*point, *line],
            Self::Symmetric {
                point_a,
                point_b,
                axis,
            } => vec![*point_a, *point_b, *axis],
        }
    }

    /// Whether this is a driving dimensional constraint (vs. geometric).
    #[must_use]
    pub fn is_driving_dimension(&self) -> bool {
        matches!(
            self,
            Self::Distance { driving: true, .. }
                | Self::Angle { driving: true, .. }
                | Self::Radius { driving: true, .. }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn horizontal_has_one_residual() {
        let c = Constraint::Horizontal(EntityId(0));
        assert_eq!(c.residual_count(), 1);
    }

    #[test]
    fn coincident_has_two_residuals() {
        let c = Constraint::Coincident(EntityId(0), EntityId(1));
        assert_eq!(c.residual_count(), 2);
    }

    #[test]
    fn distance_references_two_entities() {
        let c = Constraint::Distance {
            a: EntityId(0),
            b: EntityId(1),
            value: 5.0,
            driving: true,
        };
        assert_eq!(c.referenced_entities().len(), 2);
        assert!(c.is_driving_dimension());
    }

    #[test]
    fn geometric_constraint_is_not_driving() {
        let c = Constraint::Perpendicular(EntityId(0), EntityId(1));
        assert!(!c.is_driving_dimension());
    }
}
