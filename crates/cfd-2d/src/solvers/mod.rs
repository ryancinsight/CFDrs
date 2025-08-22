//! Numerical solvers for 2D CFD simulations
//!
//! This module contains various numerical methods for solving 2D flow problems.

pub mod fdm;
pub mod fvm;
pub mod lbm;

// Re-export main solver types
pub use fdm::{PoissonSolver, AdvectionDiffusionSolver, FdmConfig};
pub use fvm::{FvmSolver, FvmConfig, FluxScheme};
pub use lbm::{LbmSolver, LbmConfig, D2Q9};