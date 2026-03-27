//! Sketch entity types — 2D geometric primitives in sketch-local coordinates.

/// Unique identifier for a sketch entity within a single sketch.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct EntityId(pub u32);

/// A 2D point in the sketch's local coordinate system.
#[derive(Clone, Debug)]
pub struct SketchPoint {
    pub id: EntityId,
    pub x: f64,
    pub y: f64,
    pub construction: bool,
}

/// A line segment between two points.
#[derive(Clone, Debug)]
pub struct SketchLine {
    pub id: EntityId,
    pub start: EntityId,
    pub end: EntityId,
    pub construction: bool,
}

/// A circular arc defined by center, start, and end points (CCW).
#[derive(Clone, Debug)]
pub struct SketchArc {
    pub id: EntityId,
    pub center: EntityId,
    pub start: EntityId,
    pub end: EntityId,
    pub construction: bool,
}

/// A full circle defined by center point and radius.
#[derive(Clone, Debug)]
pub struct SketchCircle {
    pub id: EntityId,
    pub center: EntityId,
    pub radius: f64,
    pub construction: bool,
}

/// A NURBS spline through control points.
#[derive(Clone, Debug)]
pub struct SketchSpline {
    pub id: EntityId,
    pub control_point_ids: Vec<EntityId>,
    pub degree: usize,
    pub construction: bool,
}

/// Union of all sketch entity types.
#[derive(Clone, Debug)]
pub enum SketchEntity {
    Point(SketchPoint),
    Line(SketchLine),
    Arc(SketchArc),
    Circle(SketchCircle),
    Spline(SketchSpline),
}

impl SketchEntity {
    /// Get the entity's unique identifier.
    #[must_use]
    pub fn id(&self) -> EntityId {
        match self {
            Self::Point(e) => e.id,
            Self::Line(e) => e.id,
            Self::Arc(e) => e.id,
            Self::Circle(e) => e.id,
            Self::Spline(e) => e.id,
        }
    }

    /// Whether this entity is construction geometry (reference-only).
    #[must_use]
    pub fn is_construction(&self) -> bool {
        match self {
            Self::Point(e) => e.construction,
            Self::Line(e) => e.construction,
            Self::Arc(e) => e.construction,
            Self::Circle(e) => e.construction,
            Self::Spline(e) => e.construction,
        }
    }

    /// Get this entity as a point, if it is one.
    #[must_use]
    pub fn as_point(&self) -> Option<&SketchPoint> {
        match self {
            Self::Point(p) => Some(p),
            _ => None,
        }
    }

    /// Get this entity as a mutable point, if it is one.
    pub fn as_point_mut(&mut self) -> Option<&mut SketchPoint> {
        match self {
            Self::Point(p) => Some(p),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn entity_id_is_accessible() {
        let pt = SketchEntity::Point(SketchPoint {
            id: EntityId(42),
            x: 1.0,
            y: 2.0,
            construction: false,
        });
        assert_eq!(pt.id(), EntityId(42));
        assert!(!pt.is_construction());
    }

    #[test]
    fn construction_flag_works() {
        let line = SketchEntity::Line(SketchLine {
            id: EntityId(1),
            start: EntityId(0),
            end: EntityId(2),
            construction: true,
        });
        assert!(line.is_construction());
    }
}
