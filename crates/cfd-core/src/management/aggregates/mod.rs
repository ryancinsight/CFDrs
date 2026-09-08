//! Domain aggregates for CFD simulations.
//!
//! This module provides aggregate roots that encapsulate related entities
//! and enforce business rules and invariants.

pub mod metadata;
pub mod parameters;
pub mod problem;
pub mod simulation;
pub mod state;

pub use metadata::SimulationMetadata;
pub use parameters::PhysicalParameters;
pub use problem::ProblemAggregate;
pub use simulation::SimulationAggregate;
pub use state::SimulationState;
