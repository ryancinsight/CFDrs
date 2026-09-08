//! Pressure-velocity coupling algorithms.
//!
//! Provides configuration and result types for iterative pressure-velocity
//! coupling schemes used in incompressible flow solvers.

pub mod simple;

pub use simple::{SIMPLEConfig, SolveResult};
