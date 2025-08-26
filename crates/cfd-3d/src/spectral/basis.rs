//! Spectral basis functions and operations

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
/// Spectral basis type
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum SpectralBasis {
    /// Fourier basis for periodic domains
    Fourier,
    /// Chebyshev polynomials for non-periodic domains
    Chebyshev,
    /// Legendre polynomials
    Legendre,
}
/// Trait for basis function operations
pub trait BasisFunction<T: RealField + Copy> {
    /// Evaluate basis function at given point
    fn evaluate(&self, x: T, mode: usize) -> T;
    /// Compute derivative of basis function
    fn derivative(&self, x: T, mode: usize, order: usize) -> T;
    /// Get quadrature weights for integration
    fn quadrature_weights(&self, n: usize) -> Vec<T>;
    /// Get collocation points
    fn collocation_points(&self, n: usize) -> Vec<T>;


}
