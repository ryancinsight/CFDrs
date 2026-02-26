//! Finite Volume Method (FVM) solver for 2D CFD simulations
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

pub mod config;
pub mod flux;
pub mod geometry;
pub mod solver;

// Re-export main types
pub use config::FvmConfig;
pub use flux::{FluxScheme, FluxSchemeFactory};
pub use geometry::Face;
pub use solver::FvmSolver;
