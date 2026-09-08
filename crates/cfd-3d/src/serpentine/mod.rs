//! 3D Serpentine solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for Serpentine (sinuous) channels with:
//! - High-fidelity structured hex meshes with curvilinear mapping
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology
//! - Analysis of Dean vortices and secondary flow effects
//!
//! # Theorem — Dean Number Instability Criterion
//!
//! Secondary Dean vortices arise when
//! $De = Re\sqrt{D/(2R)}$ exceeds a critical threshold, providing the
//! governing dimensionless criterion for curved-channel transition behavior.

pub mod solver;
pub mod validation;

pub use solver::{SerpentineConfig3D, SerpentineSolution3D, SerpentineSolver3D};
pub use validation::{SerpentineValidationResult3D, SerpentineValidator3D};
