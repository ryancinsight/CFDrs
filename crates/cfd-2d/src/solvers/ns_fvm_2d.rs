//! Backward-compatibility shim for `ns_fvm_2d`.
//!
//! All types have been moved to [`super::ns_fvm`].
//! This file exists solely to avoid breaking downstream `use crate::solvers::ns_fvm_2d::*` imports.
//! **Prefer importing directly from `ns_fvm` for new code.**
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

pub use super::ns_fvm::{
    BCType, BloodModel, BoundaryCondition, FlowField2D, NavierStokesSolver2D, SIMPLEConfig,
    SolveResult, StaggeredGrid2D,
};
