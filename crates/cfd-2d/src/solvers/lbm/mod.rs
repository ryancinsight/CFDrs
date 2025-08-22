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

mod lattice;
mod collision;
mod streaming;
mod boundary;
mod solver;
mod macroscopic;

pub use lattice::{D2Q9, LatticeModel};
pub use collision::{CollisionOperator, BgkCollision};
pub use streaming::StreamingOperator;
pub use boundary::{BoundaryHandler, BoundaryType};
pub use solver::{LbmSolver, LbmConfig};
pub use macroscopic::{MacroscopicQuantities, compute_density, compute_velocity};

// Re-export for backward compatibility
pub use solver::LbmConfig as Config;