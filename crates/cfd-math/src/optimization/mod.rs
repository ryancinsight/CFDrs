//! Multi-objective optimization algorithms.
//!
//! Provides utilities for Pareto-dominance filtering and diversity-preserving
//! selection used in evolutionary multi-objective optimization (e.g., NSGA-II,
//! MOEA/D, and CFD design-space exploration).

pub mod pareto;

pub use pareto::{crowding_distances, pareto_front, ObjectiveSense};
