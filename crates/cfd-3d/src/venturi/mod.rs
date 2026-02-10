//! 3D Venturi solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for Venturi throats with:
//! - High-fidelity structured hex meshes
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology (Carreau-Yasuda, Casson)
//! - Complete validation against analytical and ISO 5167 solutions
//! - Cavitation and hemolysis prediction for biomedical applications

pub mod solver;
pub mod validation;

#[cfg(test)]
mod cavitation_hemolysis_tests;

pub use solver::{VenturiSolver3D, VenturiConfig3D, VenturiSolution3D};
pub use validation::{VenturiValidator3D, VenturiValidationResult3D};
