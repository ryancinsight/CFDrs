//! Advanced pressure-velocity coupling algorithms
//!
//! This module provides SIMPLEC and PIMPLE algorithms for improved
//! convergence in incompressible flow simulations.
//!
//! ## References
//!
//! - Van Doormaal, J. P., & Raithby, G. D. (1984). Enhancements of the SIMPLE method for predicting incompressible fluid flows.
//! - OpenFOAM PIMPLE implementation

pub mod solver;
pub mod config;

pub use config::{AlgorithmType, SimplecPimpleConfig};
pub use solver::SimplecPimpleSolver;
