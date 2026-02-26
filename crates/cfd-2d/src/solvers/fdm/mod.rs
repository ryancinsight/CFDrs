//! Finite Difference Method (FDM) solvers for 2D CFD problems.
//!
//! This module provides finite difference implementations for solving
//! various 2D fluid dynamics problems with proper separation of concerns.
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

pub mod advection_diffusion;
pub mod config;
pub mod diffusion;
pub mod linear_solver;
pub mod poisson;

// Re-export main types
pub use advection_diffusion::AdvectionDiffusionSolver;
pub use config::FdmConfig;
pub use diffusion::DiffusionSolver;
pub use linear_solver::solve_gauss_seidel;
pub use poisson::PoissonSolver;
