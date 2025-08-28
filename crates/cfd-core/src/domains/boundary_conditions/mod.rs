//! Boundary conditions domain - Physical constraints and boundary condition management.
//!
//! This module encapsulates boundary condition knowledge following DDD principles.
//! It provides abstractions for different boundary condition types and their application.

pub mod applicator;
pub mod geometry;
pub mod manager;
pub mod specification;
pub mod time_dependent;
pub mod types;

pub use applicator::BoundaryConditionApplicator;
pub use geometry::{BoundaryGeometry, BoundaryRegion};
pub use manager::BoundaryConditionManager;
pub use specification::BoundaryConditionSpec;
pub use time_dependent::{TimeDependentSpec, TimeFunctionType};
pub use types::{DirichletApplicator, NeumannApplicator, RobinApplicator};
