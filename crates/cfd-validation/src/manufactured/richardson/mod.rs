//! Richardson extrapolation and convergence analysis
//!
//! This module provides the core Richardson extrapolation algorithms for
//! grid convergence studies and error estimation in CFD validation.

pub mod core;
pub mod validation;
pub mod analysis;
pub mod types;

pub use core::*;
pub use validation::*;
pub use analysis::*;
pub use types::*;


