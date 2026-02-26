//! 2D Navier-Stokes FVM solver hierarchy.
//!
//! ## Structure
//! ```text
//! ns_fvm/
//!   mod.rs        — canonical re-exports
//!   grid.rs       — StaggeredGrid2D
//!   field.rs      — FlowField2D
//!   boundary.rs   — BCType, BoundaryCondition, BloodModel
//!   config.rs     — SIMPLEConfig, SolveResult
//!   solver.rs     — NavierStokesSolver2D + SIMPLE implementation
//! ```
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

pub mod boundary;
pub mod config;
pub mod field;
pub mod grid;
pub mod solver;

pub use boundary::{BCType, BloodModel, BoundaryCondition};
pub use config::{SIMPLEConfig, SolveResult};
pub use field::FlowField2D;
pub use grid::StaggeredGrid2D;
pub use solver::NavierStokesSolver2D;
