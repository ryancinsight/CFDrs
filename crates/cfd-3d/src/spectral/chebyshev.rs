//! Chebyshev polynomial operations for spectral methods
//!
//! Reference: Trefethen, L.N. (2000). "Spectral Methods in MATLAB"

use cfd_core::error::Result;
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;
use std::f64::consts::PI;

/// Chebyshev polynomial basis
pub struct ChebyshevPolynomial<T: RealField + Copy> {
    /// Number of collocation points
    n: usize,
    /// Collocation points (Gauss-Lobatto)
    points: Vec<T>,
    /// Differentiation matrix
    diff_matrix: DMatrix<T>,
}

impl<T: RealField + FromPrimitive + Copy> ChebyshevPolynomial<T> {
    /// Create new Chebyshev basis with n collocation points
    pub fn new(n: usize) -> Result<Self> {
        let points = Self::gauss_lobatto_points(n)?;
        let diff_matrix = Self::differentiation_matrix(n, &points)?;

        Ok(Self {
            n,
            points,
            diff_matrix,
        })
    }

    /// Get the number of collocation points
    #[must_use]
    pub fn num_points(&self) -> usize {
        self.n
    }

    /// Compute Gauss-Lobatto collocation points
    /// `x_j` = cos(πj/N) for j = 0, 1, ..., N
    fn gauss_lobatto_points(n: usize) -> Result<Vec<T>> {
        let mut points = Vec::with_capacity(n);
        let n_f64 = n as f64 - 1.0;

        for j in 0..n {
            let theta = PI * (j as f64) / n_f64;
            let x = T::from_f64(theta.cos()).ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert collocation point".into(),
                )
            })?;
            points.push(x);
        }

        Ok(points)
    }

    /// Construct Chebyshev differentiation matrix
    /// Based on Trefethen (2000), Chapter 6
    fn differentiation_matrix(n: usize, points: &[T]) -> Result<DMatrix<T>> {
        let mut d = DMatrix::zeros(n, n);
        let two = T::from_f64(2.0).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert constant".into())
        })?;

        // c_i = 2 for i = 0 or N, 1 otherwise
        let mut c = vec![T::one(); n];
        c[0] = two;
        c[n - 1] = two;

        // Fill differentiation matrix
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let sign_factor = T::from_i32(if (i + j) % 2 == 0 { 1 } else { -1 })
                        .ok_or_else(|| {
                            cfd_core::error::Error::InvalidConfiguration(
                                "Cannot convert sign factor (1 or -1) to numeric type T".into(),
                            )
                        })?;
                    let num = c[i] * sign_factor;
                    let den = c[j] * (points[i] - points[j]);
                    d[(i, j)] = num / den;
                }
            }
        }

        // Diagonal entries from sum rule
        for i in 0..n {
            let mut sum = T::zero();
            for j in 0..n {
                if i != j {
                    sum += d[(i, j)];
                }
            }
            d[(i, i)] = -sum;
        }

        Ok(d)
    }

    /// Apply differentiation operator
    #[must_use]
    pub fn differentiate(&self, u: &DVector<T>) -> DVector<T> {
        &self.diff_matrix * u
    }

    /// Get collocation points
    #[must_use]
    pub fn points(&self) -> &[T] {
        &self.points
    }

    /// Get collocation points (alias for compatibility)
    #[must_use]
    pub fn collocation_points(&self) -> &[T] {
        &self.points
    }

    /// Get differentiation matrix reference
    #[must_use]
    pub fn diff_matrix(&self) -> &DMatrix<T> {
        &self.diff_matrix
    }

    /// Get second derivative matrix (D²)
    pub fn second_derivative_matrix(&self) -> Result<DMatrix<T>> {
        // Second derivative is D * D
        Ok(&self.diff_matrix * &self.diff_matrix)
    }

    /// Interpolate function values to arbitrary point
    ///
    /// Uses barycentric Lagrange interpolation for stability
    /// Reference: Berrut & Trefethen (2004). "Barycentric Lagrange Interpolation"
    pub fn interpolate(&self, values: &[T], x: T) -> Result<T> {
        if values.len() != self.n {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "Expected {} values, got {}",
                self.n,
                values.len()
            )));
        }

        let two = T::from_f64(2.0).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert constant".into())
        })?;

        // Barycentric weights for Chebyshev points
        let mut weights = vec![T::one(); self.n];
        weights[0] = T::one() / two;
        weights[self.n - 1] = T::from_i32(if (self.n - 1).is_multiple_of(2) {
            1
        } else {
            -1
        })
        .ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert sign".into())
        })? / two;

        for (j, weight) in weights.iter_mut().enumerate().take(self.n - 1).skip(1) {
            *weight = T::from_i32(if j % 2 == 0 { 1 } else { -1 }).ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration("Cannot convert sign".into())
            })?;
        }

        // Check if x matches a grid point
        for (j, &x_j) in self.points.iter().enumerate() {
            if (x - x_j).abs() < T::from_f64(1e-14).unwrap_or_else(T::zero) {
                return Ok(values[j]);
            }
        }

        // Barycentric interpolation
        let mut numer = T::zero();
        let mut denom = T::zero();

        for j in 0..self.n {
            let term = weights[j] / (x - self.points[j]);
            numer += term * values[j];
            denom += term;
        }

        Ok(numer / denom)
    }
}

// ChebyshevDifferentiation struct removed as it was redundant.
// Users can call differentiate() directly on ChebyshevPolynomial instances.

#[cfg(test)]
#[path = "chebyshev_tests.rs"]
mod chebyshev_tests;
