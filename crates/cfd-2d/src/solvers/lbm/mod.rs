//! Lattice Boltzmann Method (LBM) solvers for 2D fluid dynamics.
//!
//! This module provides LBM implementations for solving incompressible
//! Navier-Stokes equations using the collision-streaming approach.
//!
//! Features:
//! - D2Q9 lattice model (2D, 9 velocities)
//! - BGK collision operator
//! - Various boundary conditions (bounce-back, velocity, pressure)
//! - Parallel processing support

mod boundary;
mod collision;
mod lattice;
mod macroscopic;
mod solver;
mod streaming;

pub use boundary::{BoundaryHandler, BoundaryType};
pub use collision::{BgkCollision, CollisionOperator};
pub use lattice::{LatticeModel, D2Q9};
pub use macroscopic::{compute_density, compute_velocity, MacroscopicQuantities};
pub use solver::{LbmConfig, LbmSolver};
pub use streaming::StreamingOperator;

// Re-export for backward compatibility
pub use solver::LbmConfig as Config;
