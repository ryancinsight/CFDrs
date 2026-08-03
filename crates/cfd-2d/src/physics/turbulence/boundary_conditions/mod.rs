//! Turbulence boundary condition implementations
//!
//! This module provides physically accurate boundary conditions for turbulence models:
//! - Wall boundary conditions for k, ε, ω, ν̃
//! - Inlet/outlet boundary conditions for turbulence quantities
//! - Wall distance calculation for wall functions
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

mod manager;
mod wall;

#[cfg(test)]
mod tests;

pub use manager::{TurbulenceBoundaryCondition, TurbulenceBoundaryManager};
