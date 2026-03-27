//! Sketch undo/redo commands — entity and constraint add/remove operations.

use crate::domain::sketch::constraint::{Constraint, ConstraintId};
use crate::domain::sketch::entity::{EntityId, SketchEntity};
use crate::domain::sketch::sketch::SketchHandle;

/// Command to add a sketch entity.
#[derive(Clone, Debug)]
pub struct AddSketchEntityCommand {
    /// Which sketch to modify.
    pub sketch_handle: SketchHandle,
    /// The entity to add.
    pub entity: SketchEntity,
}

/// Command to remove a sketch entity.
#[derive(Clone, Debug)]
pub struct RemoveSketchEntityCommand {
    /// Which sketch to modify.
    pub sketch_handle: SketchHandle,
    /// ID of the entity to remove.
    pub entity_id: EntityId,
    /// The removed entity (for undo).
    pub removed: Option<SketchEntity>,
}

/// Command to add a constraint.
#[derive(Clone, Debug)]
pub struct AddConstraintCommand {
    /// Which sketch to modify.
    pub sketch_handle: SketchHandle,
    /// The constraint to add.
    pub constraint: Constraint,
    /// The assigned constraint ID (for undo).
    pub assigned_id: Option<ConstraintId>,
}

/// Command to remove a constraint.
#[derive(Clone, Debug)]
pub struct RemoveConstraintCommand {
    /// Which sketch to modify.
    pub sketch_handle: SketchHandle,
    /// ID of the constraint to remove.
    pub constraint_id: ConstraintId,
    /// The removed constraint (for undo).
    pub removed: Option<Constraint>,
}
