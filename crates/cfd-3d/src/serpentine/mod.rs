//! 3D Serpentine solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for Serpentine (sinuous) channels with:
//! - High-fidelity structured hex meshes with curvilinear mapping
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology
//! - Analysis of Dean vortices and secondary flow effects

pub mod solver;
pub mod validation;

pub use solver::{SerpentineSolver3D, SerpentineConfig3D, SerpentineSolution3D};
pub use validation::{SerpentineValidator3D, SerpentineValidationResult3D};
