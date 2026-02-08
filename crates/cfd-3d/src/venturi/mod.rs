//! 3D Venturi solvers with FEM validation
//!
//! This module implements full 3D CFD simulations for Venturi throats with:
//! - High-fidelity structured hex meshes
//! - Finite Element Method (FEM) for Navier-Stokes equations
//! - Blood flow with non-Newtonian rheology (Carreau-Yasuda, Casson)
//! - Complete validation against analytical and ISO 5167 solutions

pub mod solver;
pub use solver::{VenturiSolver3D, VenturiConfig3D, VenturiSolution3D};
