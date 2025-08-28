//! Collision operators for Lattice Boltzmann Method
//!
//! Clean separation of collision models following SOLID principles

mod bgk;
mod forcing;
mod mrt;
mod regularized;
mod traits;

pub use bgk::BgkCollision;
pub use traits::CollisionOperator;
