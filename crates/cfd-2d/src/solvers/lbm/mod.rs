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
//!
//! # Theorem
//! The Lattice Boltzmann Method (LBM) recovers the macroscopic Navier-Stokes equations
//! in the low Mach number limit.
//!
//! **Proof sketch**:
//! Through the Chapman-Enskog expansion, the discrete Boltzmann equation with the BGK
//! collision operator can be expanded in powers of the Knudsen number ($Kn$).
//! At $O(Kn^0)$, the Euler equations are recovered. At $O(Kn^1)$, the viscous stress
//! tensor emerges, yielding the weakly compressible Navier-Stokes equations.
//! The kinematic viscosity is related to the relaxation time $\tau$ by $\nu = c_s^2 (\tau - 0.5)\Delta t$.

mod boundary;
mod collision;
mod lattice;
mod macroscopic;
/// Multiphase LBM models (Shan-Chen pseudopotential).
pub mod multiphase;
mod scalar_boundary;
mod solver;
mod streaming;

pub use boundary::{BoundaryHandler, BoundaryType};
pub use collision::{
    BgkCollision, CarreauYasudaBgk, CollisionOperator, MrtCollision, RelaxationMatrix,
};
pub use lattice::{LatticeModel, D2Q9};
pub use macroscopic::{
    compute_density, compute_kinetic_energy, compute_stress_tensor, compute_velocity,
    compute_vorticity, MacroscopicQuantities,
};
pub use multiphase::ShanChenMultiphase;
pub use solver::{LbmConfig, LbmSolver};
pub use streaming::StreamingOperator;
