//! Numerical schemes for 2D CFD simulations
//!
//! This module provides various discretization schemes for spatial and temporal
//! derivatives in 2D computational fluid dynamics.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

// Re-export submodules
pub mod upwind;
pub mod central;
pub mod tvd;
pub mod weno;
pub mod time_integration;
pub mod grid;
pub mod constants;

// Re-export main types
pub use upwind::{FirstOrderUpwind, SecondOrderUpwind};
pub use central::{CentralDifference, FourthOrderCentral};
pub use tvd::{MUSCLScheme, QUICKScheme, FluxLimiter};
pub use weno::WENO5;
pub use time_integration::{TimeScheme, TimeIntegrator};
pub use grid::Grid2D;

/// Spatial discretization scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SpatialScheme {
    /// First-order upwind
    FirstOrderUpwind,
    /// Second-order central differencing
    CentralDifference,
    /// Second-order upwind
    SecondOrderUpwind,
    /// QUICK (Quadratic Upstream Interpolation)
    Quick,
    /// Third-order MUSCL
    Muscl,
    /// Fifth-order WENO
    Weno5,
    /// Fourth-order explicit central difference
    FourthOrderCentral,
}

/// Trait for spatial discretization schemes
pub trait SpatialDiscretization<T: RealField> {
    /// Compute spatial derivative
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T;
    
    /// Get scheme order of accuracy
    fn order(&self) -> usize;
    
    /// Check if scheme is conservative
    fn is_conservative(&self) -> bool;
}