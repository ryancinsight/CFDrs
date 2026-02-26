//! Spectral basis functions and operations
//!
//! # Theorem — Fourier Spectral Convergence (Canuto et al. 2006)
//!
//! For a periodic function $u$ with Fourier coefficients $\hat{u}_k$:
//! - If $u \in C^m$, then $|\hat{u}_k| = O(k^{-m})$ (algebraic decay).
//! - If $u$ is analytic in a strip of width $\sigma$, then $|\hat{u}_k| = O(e^{-\sigma|k|})$
//!   (exponential/spectral convergence).
//!
//! # Theorem — Chebyshev Minimax Error Bound (Trefethen 2013)
//!
//! The Chebyshev interpolant $p_N$ of degree $N$ satisfies
//!
//! ```text
//! ‖u − p_N‖_∞ ≤ (1 + Λ_N) E_N(u)
//! ```
//!
//! where $E_N(u)$ is the minimax polynomial error and the Lebesgue
//! constant $\Lambda_N = O(\log N)$ for Chebyshev–Gauss–Lobatto points.
//! For analytic $u$, $E_N(u) = O(e^{-cN})$, so convergence is spectral.

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
pub trait BasisFunction<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Evaluate basis function at given point
    fn evaluate(&self, x: T, mode: usize) -> T;

    /// Compute derivative of basis function
    fn derivative(&self, x: T, mode: usize, order: usize) -> T;

    /// Get quadrature weights for integration
    fn quadrature_weights(&self, n: usize) -> Vec<T>;

    /// Get collocation points
    fn collocation_points(&self, n: usize) -> Vec<T>;
}
