//! Collision operators for Lattice Boltzmann Method
//!
//! Clean separation of collision models following SOLID principles
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

mod bgk;
mod mrt;
mod regularized;
mod traits;

pub use bgk::BgkCollision;
pub use mrt::{MrtCollision, RelaxationMatrix};
pub use traits::CollisionOperator;
