//! Richardson extrapolation and convergence analysis
//!
//! This module provides the core Richardson extrapolation algorithms for
//! grid convergence studies and error estimation in CFD validation.

pub mod analysis;
pub mod core;
pub mod types;
pub mod validation;

pub use analysis::*;
pub use core::*;
pub use types::*;
pub use validation::*;

