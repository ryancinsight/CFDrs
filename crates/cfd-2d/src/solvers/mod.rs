//! Numerical solvers for 2D CFD simulations
//!
//! This module contains various numerical methods for solving 2D flow problems.

pub mod accelerated;
pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod ns_fvm_2d;
pub mod poiseuille;
pub mod serpentine_flow;
pub mod simd_kernels;
pub mod simple;
pub mod venturi_flow;

// Re-export main solver types
pub use fdm::{AdvectionDiffusionSolver, DiffusionSolver, FdmConfig, PoissonSolver};
pub use fvm::{FluxScheme, FvmConfig, FvmSolver};
pub use lbm::{LbmConfig, LbmSolver, D2Q9};
pub use ns_fvm_2d::NavierStokesSolver2D;
pub use poiseuille::{BloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
pub use simple::SimpleAlgorithm;
