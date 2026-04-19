//! Numerical schemes for 2D CFD simulations
//!
//! This module provides various discretization schemes for spatial and temporal
//! derivatives in 2D computational fluid dynamics.
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

// Re-export submodules
pub mod central;
pub mod constants;
pub mod grid;
pub mod time;
pub mod tvd;
#[cfg(test)]
/// TVD scheme regression and theorem-validation tests.
pub mod tvd_tests;
pub mod upwind;
pub mod weno;
pub mod weno_constants;
pub(crate) mod weno_helpers;
#[cfg(test)]
pub mod weno_tests;
pub mod weno_z;

// Re-export main types
pub use central::{CentralDifference, FourthOrderCentral};
pub use grid::Grid2D;
pub use time::{TimeIntegrator, TimeScheme};
pub use tvd::{FluxLimiter, MUSCLScheme, QUICKScheme};
pub use upwind::{FirstOrderUpwind, SecondOrderUpwind};
pub use weno::{WENO5, WENO9};
pub use weno_z::WENOZ5;

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
    /// Fifth-order WENO-Z
    WenoZ5,
    /// Ninth-order WENO
    Weno9,
    /// Fourth-order explicit central difference
    FourthOrderCentral,
}

/// Trait for spatial discretization schemes
pub trait SpatialDiscretization<T: RealField + Copy> {
    /// Compute the spatial derivative on a cell-centered grid.
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
