//! Geometric primitives, exact predicates, and mesh builders.

pub mod plane;
pub mod aabb;
pub mod normal;
pub mod measure;
pub mod distance;
pub mod predicates;

pub use plane::Plane;
pub use aabb::Aabb;
pub use predicates::Orientation;

pub mod venturi;
pub mod serpentine;
pub mod branching;
