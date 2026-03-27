//! Sketch struct — the core 2D parametric sketch data model.

use std::ops::Range;

use super::constraint::{Constraint, ConstraintId};
use super::entity::{EntityId, SketchEntity, SketchPoint};
use super::work_plane::WorkPlane;

/// Unique identifier for a sketch in the project.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SketchId(pub u32);

/// Handle into the project document's sketch arena.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SketchHandle(pub usize);

/// A 2D parametric sketch living on a work plane.
#[derive(Clone, Debug)]
pub struct Sketch {
    pub id: SketchId,
    pub name: String,
    pub work_plane: WorkPlane,
    entities: Vec<SketchEntity>,
    constraints: Vec<(ConstraintId, Constraint)>,
    next_entity_id: u32,
    next_constraint_id: u32,
}

impl Sketch {
    /// Create a new empty sketch on the given work plane.
    #[must_use]
    pub fn new(id: SketchId, name: String, work_plane: WorkPlane) -> Self {
        Self {
            id,
            name,
            work_plane,
            entities: Vec::new(),
            constraints: Vec::new(),
            next_entity_id: 0,
            next_constraint_id: 0,
        }
    }

    /// Add an entity to the sketch.
    pub fn add_entity(&mut self, entity: SketchEntity) {
        self.entities.push(entity);
    }

    /// Allocate the next available entity ID.
    pub fn next_entity_id(&mut self) -> EntityId {
        let id = EntityId(self.next_entity_id);
        self.next_entity_id += 1;
        id
    }

    /// Remove an entity by ID.
    pub fn remove_entity(&mut self, id: EntityId) {
        self.entities.retain(|e| e.id() != id);
    }

    /// Look up an entity by ID.
    #[must_use]
    pub fn entity(&self, id: EntityId) -> Option<&SketchEntity> {
        self.entities.iter().find(|e| e.id() == id)
    }

    /// Mutable access to an entity by ID.
    pub fn entity_mut(&mut self, id: EntityId) -> Option<&mut SketchEntity> {
        self.entities.iter_mut().find(|e| e.id() == id)
    }

    /// All entities in the sketch.
    #[must_use]
    pub fn entities(&self) -> &[SketchEntity] {
        &self.entities
    }

    /// Look up a point entity by ID.
    #[must_use]
    pub fn point(&self, id: EntityId) -> Option<&SketchPoint> {
        self.entity(id).and_then(|e| e.as_point())
    }

    /// Mutable access to a point entity.
    pub fn point_mut(&mut self, id: EntityId) -> Option<&mut SketchPoint> {
        self.entity_mut(id).and_then(|e| e.as_point_mut())
    }

    /// Add a constraint. Returns the assigned ID.
    pub fn add_constraint(&mut self, constraint: Constraint) -> ConstraintId {
        let id = ConstraintId(self.next_constraint_id);
        self.next_constraint_id += 1;
        self.constraints.push((id, constraint));
        id
    }

    /// Remove a constraint by ID.
    pub fn remove_constraint(&mut self, id: ConstraintId) {
        self.constraints.retain(|(cid, _)| *cid != id);
    }

    /// All constraints in the sketch.
    #[must_use]
    pub fn constraints(&self) -> &[(ConstraintId, Constraint)] {
        &self.constraints
    }

    /// Collect all free DOF parameters into a flat vector.
    ///
    /// For each point: `[x, y]`. For each circle: appends `[radius]`.
    /// The parameter map provides the index range for each entity.
    #[must_use]
    pub fn parameter_vector(&self) -> Vec<f64> {
        let mut params = Vec::new();
        for entity in &self.entities {
            match entity {
                SketchEntity::Point(p) => {
                    params.push(p.x);
                    params.push(p.y);
                }
                SketchEntity::Circle(c) => {
                    params.push(c.radius);
                }
                _ => {}
            }
        }
        params
    }

    /// Write solved parameters back into entities.
    pub fn apply_parameters(&mut self, params: &[f64]) {
        let mut idx = 0;
        for entity in &mut self.entities {
            match entity {
                SketchEntity::Point(p) if idx + 1 < params.len() => {
                    p.x = params[idx];
                    p.y = params[idx + 1];
                    idx += 2;
                }
                SketchEntity::Circle(c) if idx < params.len() => {
                    c.radius = params[idx];
                    idx += 1;
                }
                _ => {}
            }
        }
    }

    /// Map each entity to its parameter index range.
    #[must_use]
    pub fn parameter_map(&self) -> Vec<(EntityId, Range<usize>)> {
        let mut map = Vec::new();
        let mut idx = 0;
        for entity in &self.entities {
            match entity {
                SketchEntity::Point(p) => {
                    map.push((p.id, idx..idx + 2));
                    idx += 2;
                }
                SketchEntity::Circle(c) => {
                    map.push((c.id, idx..idx + 1));
                    idx += 1;
                }
                _ => {}
            }
        }
        map
    }

    /// Total number of entities.
    #[must_use]
    pub fn entity_count(&self) -> usize {
        self.entities.len()
    }

    /// Total number of constraints.
    #[must_use]
    pub fn constraint_count(&self) -> usize {
        self.constraints.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::entity::SketchLine;

    fn test_sketch() -> Sketch {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let ln = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p0, x: 0.0, y: 0.0, construction: false,
        }));
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id: p1, x: 3.0, y: 4.0, construction: false,
        }));
        sk.add_entity(SketchEntity::Line(SketchLine {
            id: ln, start: p0, end: p1, construction: false,
        }));
        sk
    }

    #[test]
    fn parameter_vector_has_correct_length() {
        let sk = test_sketch();
        let params = sk.parameter_vector();
        assert_eq!(params.len(), 4); // 2 points × 2 DOFs
    }

    #[test]
    fn apply_parameters_updates_points() {
        let mut sk = test_sketch();
        sk.apply_parameters(&[10.0, 20.0, 30.0, 40.0]);
        let p0 = sk.point(EntityId(0)).expect("exists");
        assert_eq!(p0.x, 10.0);
        assert_eq!(p0.y, 20.0);
    }

    #[test]
    fn parameter_map_ranges() {
        let sk = test_sketch();
        let map = sk.parameter_map();
        assert_eq!(map.len(), 2); // 2 points (line has no params)
        assert_eq!(map[0].1, 0..2);
        assert_eq!(map[1].1, 2..4);
    }

    #[test]
    fn add_and_remove_constraint() {
        let mut sk = test_sketch();
        let cid = sk.add_constraint(Constraint::Horizontal(EntityId(2)));
        assert_eq!(sk.constraint_count(), 1);
        sk.remove_constraint(cid);
        assert_eq!(sk.constraint_count(), 0);
    }
}
