//! Numerical solvers for 2D CFD simulations
//!
//! This module contains various numerical methods for solving 2D flow problems.

pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod simd_kernels;

// Re-export main solver types
pub use fdm::{AdvectionDiffusionSolver, FdmConfig, PoissonSolver};
pub use fvm::{FluxScheme, FvmConfig, FvmSolver};
pub use lbm::{LbmConfig, LbmSolver, D2Q9};
