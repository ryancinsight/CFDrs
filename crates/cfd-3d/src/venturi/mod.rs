//! 3D Venturi solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for Venturi throats with:
//! - High-fidelity structured hex meshes
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology (Carreau-Yasuda, Casson)
//! - Complete validation against analytical and ISO 5167 solutions
//! - Cavitation and hemolysis prediction for biomedical applications
//!
//! # Theorem — Venturi Energy Relation
//!
//! For incompressible steady flow, continuity and Bernoulli imply a throat
//! pressure drop linked to area contraction; in viscous realizations,
//! measured drop is bounded below by the inviscid Bernoulli prediction.

pub mod solver;
pub mod validation;

#[cfg(test)]
mod cavitation_hemolysis_tests;

pub use solver::{VenturiConfig3D, VenturiSolution3D, VenturiSolver3D};
pub use validation::{VenturiValidationResult3D, VenturiValidator3D};
