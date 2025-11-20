//! Numerical solvers for 2D CFD simulations
//!
//! This module contains various numerical methods for solving 2D flow problems.

pub mod accelerated;
pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod simd_kernels;
pub mod simple;

// Re-export main solver types
pub use fdm::{AdvectionDiffusionSolver, DiffusionSolver, FdmConfig, PoissonSolver};
pub use fvm::{FluxScheme, FvmConfig, FvmSolver};
pub use lbm::{LbmConfig, LbmSolver, D2Q9};
pub use simple::SimpleAlgorithm;
