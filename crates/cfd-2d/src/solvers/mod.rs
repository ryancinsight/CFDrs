//! Numerical solvers for 2D CFD simulations
//!
//! This module contains various numerical methods for solving 2D flow problems.

pub mod accelerated;
pub mod cavity_solver;
pub mod fdm;
pub mod fvm;
pub mod lbm;
pub mod ns_fvm_2d;
pub mod scalar_transport_2d;
pub mod bifurcation_flow;
pub mod poiseuille;
pub mod serpentine_flow;
pub mod simd_kernels;
pub mod simple;
pub mod venturi_flow;

// Re-export main solver types
pub use bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
pub use fdm::{AdvectionDiffusionSolver, DiffusionSolver, FdmConfig, PoissonSolver};
pub use fvm::{FluxScheme, FvmConfig, FvmSolver};
pub use lbm::{LbmConfig, LbmSolver, D2Q9};
pub use poiseuille::{BloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
pub use simple::SimpleAlgorithm;
