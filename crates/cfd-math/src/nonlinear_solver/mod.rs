//! # Non-linear Solvers
//!
//! Algorithms for finding roots and fixed-points of non-linear systems:
//!
//! - [`anderson`]: Anderson Acceleration (Type-II MGS-QR variant, VecDeque history)
//! - [`jfnk`]: Jacobian-Free Newton-Krylov with finite-difference JvP and GMRES inner solver

pub mod anderson;
pub mod jfnk;

pub use anderson::{AndersonAccelerator, AndersonConfig, AndersonMethod};
pub use jfnk::{JfnkConfig, JfnkConvergence, JfnkSolver};
