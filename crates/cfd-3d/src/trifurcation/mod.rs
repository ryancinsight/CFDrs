//! 3D Trifurcation flow modeling and simulation
//!
//! This module provides tools for simulating flow in vessels branching into
//! three daughter vessels, common in microfluidic junctions and specialized
//! vascular regions.

pub mod geometry;
pub mod solver;

pub use geometry::TrifurcationGeometry3D;
pub use solver::{TrifurcationConfig3D, TrifurcationSolution3D, TrifurcationSolver3D};

/// Mesh for trifurcation (reuses bifurcation mesh structure)
pub type TrifurcationMesh<T> = crate::bifurcation::BifurcationMesh<T>;

#[cfg(test)]
mod tests;
