//! Level set method for interface tracking
//!
//! This module implements the level set method for tracking interfaces
//! between two fluids or phases.

pub mod config;
pub mod solver;
// Re-export main types
pub use config::{LevelSetConfig, EPSILON_DIVISION, HALF_VALUE, POWER_EXPONENT_TWO};
pub use solver::LevelSetSolver;
