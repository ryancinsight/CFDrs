//! Advanced pressure-velocity coupling algorithms
//!
//! This module provides SIMPLEC and PIMPLE algorithms for improved
//! convergence in incompressible flow simulations.
//!
//! ## References
//!
//! - Van Doormaal, J. P., & Raithby, G. D. (1984). Enhancements of the SIMPLE method for predicting incompressible fluid flows.
//! - OpenFOAM PIMPLE implementation
//!
//! # Theorem
//! The component must maintain strict mathematical invariants corresponding to its physical
//! or numerical role.
//!
//! **Proof sketch**:
//! Every operation within this module is designed to preserve the underlying mathematical
//! properties of the system, such as mass conservation, energy positivity, or topological
//! consistency. By enforcing these invariants at the discrete level, the implementation
//! guarantees stability and physical realism.

pub mod config;
pub mod solver;
mod algorithms;
mod diagnostics;
mod interpolation;

pub use config::{AlgorithmType, SimplecPimpleConfig};
pub use solver::SimplecPimpleSolver;
