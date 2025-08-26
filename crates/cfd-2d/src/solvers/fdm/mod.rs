//! Finite Difference Method (FDM) solvers for 2D CFD problems.
//!
//! This module provides finite difference implementations for solving
//! various 2D fluid dynamics problems with proper separation of concerns.

pub mod advection_diffusion;
pub mod config;
pub mod linear_solver;
pub mod poisson;
// Re-export main types
pub use advection_diffusion::AdvectionDiffusionSolver;
pub use config::FdmConfig;
pub use linear_solver::solve_gauss_seidel;
pub use poisson::PoissonSolver;
