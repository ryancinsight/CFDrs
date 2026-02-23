//! Numerical schemes for 2D CFD simulations
//!
//! This module provides various discretization schemes for spatial and temporal
//! derivatives in 2D computational fluid dynamics.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

// Re-export submodules
pub mod central;
pub mod constants;
pub mod grid;
pub mod time;
pub mod tvd;
pub mod upwind;
pub mod weno;
pub mod weno_constants;
#[cfg(test)]
pub mod tvd_tests;
#[cfg(test)]
pub mod weno_tests;

// Re-export main types
pub use central::{CentralDifference, FourthOrderCentral};
pub use grid::Grid2D;
pub use time::{TimeIntegrator, TimeScheme};
pub use tvd::{FluxLimiter, MUSCLScheme, QUICKScheme};
pub use upwind::{FirstOrderUpwind, SecondOrderUpwind};
pub use weno::{WENO5, WENO9};

/// Spatial discretization scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SpatialScheme {
    /// First-order upwind
    FirstOrderUpwind,
    /// Second-order central differencing
    CentralDifference,
    /// Second-order upwind
    SecondOrderUpwind,
    /// Quadratic Upstream Interpolation for Convective Kinematics
    QuadraticUpstreamInterpolation,
    /// Third-order MUSCL
    Muscl,
    /// Fifth-order WENO
    Weno5,
    /// Ninth-order WENO
    Weno9,
    /// Fourth-order explicit central difference
    FourthOrderCentral,
}

/// Trait for spatial discretization schemes
pub trait SpatialDiscretization<T: RealField + Copy> {
    /// Compute spatial derivative (for backward compatibility)
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T;

    /// Get scheme order of accuracy
    fn order(&self) -> usize;

    /// Check if scheme is conservative
    fn is_conservative(&self) -> bool;

    /// Compute CFL condition for advection: c*dt/dx
    /// Returns the maximum CFL number for stability
    fn cfl_limit(&self) -> f64;

    /// Check if scheme is stable for given CFL number
    fn is_stable(&self, cfl: f64) -> bool {
        cfl <= self.cfl_limit()
    }

    /// Compute von Neumann amplification factor for given wavenumber and CFL
    /// For 1D advection: G(k) = 1 - i * CFL * sin(k*dx)
    fn amplification_factor(&self, k: f64, cfl: f64) -> num_complex::Complex<f64> {
        let sin_kdx = (k * std::f64::consts::PI).sin(); // Assuming dx=1 for normalized analysis
        num_complex::Complex::new(1.0, 0.0) - num_complex::Complex::new(0.0, cfl * sin_kdx)
    }
}

/// Trait for face reconstruction schemes used in finite volume methods
pub trait FaceReconstruction<T: RealField + Copy> {
    /// Reconstruct scalar value at x-face (between cells i and i+1)
    fn reconstruct_face_value_x(
        &self,
        phi: &Grid2D<T>,     // The scalar field being transported
        velocity_at_face: T, // The velocity normal to the face
        i: usize,            // Index of the "left" cell
        j: usize,
    ) -> T;

    /// Reconstruct scalar value at y-face (between cells j and j+1)
    fn reconstruct_face_value_y(
        &self,
        phi: &Grid2D<T>,     // The scalar field being transported
        velocity_at_face: T, // The velocity normal to the face
        i: usize,
        j: usize, // Index of the "bottom" cell
    ) -> T;

    /// Get scheme order of accuracy
    fn order(&self) -> usize;
}
