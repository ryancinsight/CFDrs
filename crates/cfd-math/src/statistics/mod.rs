//! Statistical utilities for design-space analysis.
//!
//! Currently provides Pareto-dominance filtering and NSGA-II crowding distance
//! for multi-objective optimisation post-processing.

pub mod pareto;

pub use pareto::{crowding_distances, pareto_front_nd};
