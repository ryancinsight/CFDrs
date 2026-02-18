//! Geometric primitives and mesh builders.

pub mod plane;
pub mod aabb;
pub mod normal;
pub mod measure;
pub mod distance;

pub use plane::Plane;
pub use aabb::Aabb;

pub mod venturi;
pub mod serpentine;
pub mod branching;
