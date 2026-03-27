//! Sketch tool service — creates sketch geometry from user interactions.

use crate::domain::sketch::entity::{EntityId, SketchCircle, SketchEntity, SketchLine, SketchPoint};
use crate::domain::sketch::sketch::Sketch;

/// Service for creating sketch geometry from tool interactions.
pub struct SketchToolService;

impl SketchToolService {
    /// Create a line from two click positions.
    ///
    /// Creates two points and a line entity. Returns `[p0_id, p1_id, line_id]`.
    pub fn create_line(sketch: &mut Sketch, p0: (f64, f64), p1: (f64, f64)) -> Vec<EntityId> {
        let id0 = sketch.next_entity_id();
        let id1 = sketch.next_entity_id();
        let lid = sketch.next_entity_id();

        sketch.add_entity(SketchEntity::Point(SketchPoint {
            id: id0, x: p0.0, y: p0.1, construction: false,
        }));
        sketch.add_entity(SketchEntity::Point(SketchPoint {
            id: id1, x: p1.0, y: p1.1, construction: false,
        }));
        sketch.add_entity(SketchEntity::Line(SketchLine {
            id: lid, start: id0, end: id1, construction: false,
        }));

        vec![id0, id1, lid]
    }

    /// Create a circle from center and a point on the circumference.
    ///
    /// Returns `[center_id, circle_id]`.
    pub fn create_circle(
        sketch: &mut Sketch,
        center: (f64, f64),
        radius_pt: (f64, f64),
    ) -> Vec<EntityId> {
        let cid = sketch.next_entity_id();
        let circle_id = sketch.next_entity_id();

        let dx = radius_pt.0 - center.0;
        let dy = radius_pt.1 - center.1;
        let radius = (dx * dx + dy * dy).sqrt();

        sketch.add_entity(SketchEntity::Point(SketchPoint {
            id: cid, x: center.0, y: center.1, construction: false,
        }));
        sketch.add_entity(SketchEntity::Circle(SketchCircle {
            id: circle_id, center: cid, radius, construction: false,
        }));

        vec![cid, circle_id]
    }

    /// Auto-detect nearby points and add coincident constraints.
    ///
    /// Returns the IDs of any new constraints added.
    pub fn auto_constrain_coincident(
        sketch: &mut Sketch,
        point_id: EntityId,
        snap_radius: f64,
    ) -> Vec<crate::domain::sketch::constraint::ConstraintId> {
        let Some(pt) = sketch.point(point_id) else { return vec![] };
        let px = pt.x;
        let py = pt.y;

        let nearby: Vec<EntityId> = sketch
            .entities()
            .iter()
            .filter_map(|e| {
                if let SketchEntity::Point(other) = e {
                    if other.id != point_id {
                        let dx = other.x - px;
                        let dy = other.y - py;
                        if (dx * dx + dy * dy).sqrt() < snap_radius {
                            return Some(other.id);
                        }
                    }
                }
                None
            })
            .collect();

        let mut constraints = Vec::new();
        for other_id in nearby {
            let cid = sketch.add_constraint(
                crate::domain::sketch::constraint::Constraint::Coincident(point_id, other_id),
            );
            constraints.push(cid);
        }
        constraints
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;

    #[test]
    fn create_line_adds_three_entities() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let ids = SketchToolService::create_line(&mut sk, (0.0, 0.0), (1.0, 1.0));
        assert_eq!(ids.len(), 3);
        assert_eq!(sk.entity_count(), 3);
    }

    #[test]
    fn create_circle_computes_radius() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let ids = SketchToolService::create_circle(&mut sk, (0.0, 0.0), (3.0, 4.0));
        assert_eq!(ids.len(), 2);
        if let Some(SketchEntity::Circle(c)) = sk.entity(ids[1]) {
            assert!((c.radius - 5.0).abs() < 1e-12);
        } else {
            panic!("expected circle");
        }
    }
}
