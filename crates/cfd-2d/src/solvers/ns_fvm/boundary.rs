//! Boundary condition types and blood rheology model selector.
//!
//! All types are re-exported from their canonical homes in `cfd-core` and `cfd-math`.
//! This module is a thin pass-through; no types are defined here.
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

// Canonical boundary condition types from cfd-core
pub use cfd_core::physics::boundary::{BoundaryCondition, FundamentalBCType as BCType};

// Canonical blood rheology dispatch from cfd-core
pub use cfd_core::physics::fluid::BloodModel;
