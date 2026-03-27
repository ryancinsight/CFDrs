//! Sketch domain module — 2D parametric sketch entities, constraints, and work planes.

pub mod constraint;
pub mod dof;
pub mod entity;
#[allow(clippy::module_inception)]
pub mod sketch;
pub mod work_plane;

pub use constraint::{Constraint, ConstraintId};
pub use dof::{DofAnalysis, DofStatus};
pub use entity::{EntityId, SketchEntity};
pub use sketch::{Sketch, SketchHandle, SketchId};
pub use work_plane::WorkPlane;
