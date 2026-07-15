//! Chebyshev polynomial operations for spectral methods
//!
//! # Theorem — Gauss–Lobatto Exactness (Bernardi & Maday 1997)
//!
//! The $(N+1)$-point Chebyshev–Gauss–Lobatto quadrature
//!
//! ```text
//! x_j = cos(π j / N),   j = 0, 1, …, N
//! ```
//!
//! integrates polynomials of degree $\leq 2N - 1$ exactly against the
//! Chebyshev weight $w(x) = (1 - x^2)^{-1/2}$.
//!
//! # Theorem — Spectral Differentiation Convergence (Trefethen 2000)
//!
//! For an analytic function $u$, the Chebyshev differentiation matrix $D_N$
//! satisfies
//!
//! ```text
//! ‖u' − D_N u‖_∞ = O(e^{−cN})
//! ```
//!
//! where $c > 0$ depends on the analyticity strip width of $u$.
//!
//! Reference: Trefethen, L.N. (2000). "Spectral Methods in MATLAB"

use crate::scalar;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};
use leto::{Array1, Array2};
use std::f64::consts::PI;

/// Chebyshev polynomial basis
pub struct ChebyshevPolynomial<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Number of collocation points
    n: usize,
    /// Collocation points (Gauss-Lobatto)
    points: Vec<T>,
    /// Differentiation matrix
    diff_matrix: Array2<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> ChebyshevPolynomial<T> {
    /// Create new Chebyshev basis with n collocation points
    pub fn new(n: usize) -> Result<Self> {
        if n < 2 {
            return Err(cfd_core::error::Error::InvalidConfiguration(format!(
                "ChebyshevPolynomial requires at least 2 collocation points, got {n}"
            )));
        }

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
            let x = scalar::from_f64::<T>(theta.cos());
            points.push(x);
        }

        Ok(points)
    }

    /// Construct Chebyshev differentiation matrix
    /// Based on Trefethen (2000), Chapter 6
    fn differentiation_matrix(n: usize, points: &[T]) -> Result<Array2<T>> {
        let mut d = Array2::from_elem([n, n], scalar::zero::<T>());
        let two = scalar::from_f64::<T>(2.0);

        // c_i = 2 for i = 0 or N, 1 otherwise
        let mut c = vec![scalar::one::<T>(); n];
        c[0] = two;
        c[n - 1] = two;

        // Fill differentiation matrix
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let sign_factor =
                        scalar::from_f64::<T>(if (i + j) % 2 == 0 { 1.0 } else { -1.0 });
                    let num = c[i] * sign_factor;
                    let den = c[j] * (points[i] - points[j]);
                    d[[i, j]] = num / den;
                }
            }
        }

        // Diagonal entries from sum rule
        for i in 0..n {
            let mut sum = scalar::zero::<T>();
            for j in 0..n {
                if i != j {
                    sum += d[[i, j]];
                }
            }
            d[[i, i]] = -sum;
        }

        Ok(d)
    }

    /// Apply differentiation operator
    ///
    /// # Errors
    ///
    /// Returns [`cfd_core::error::Error::DimensionMismatch`] when `u` does not
    /// contain exactly one value per collocation point.
    pub fn differentiate(&self, u: &Array1<T>) -> Result<Array1<T>> {
        Self::mat_vec(&self.diff_matrix, u)
    }

    /// Get collocation points
    #[must_use]
    pub fn points(&self) -> &[T] {
        &self.points
    }

    /// Get differentiation matrix reference
    #[must_use]
    pub fn diff_matrix(&self) -> &Array2<T> {
        &self.diff_matrix
    }

    /// Get second derivative matrix (D²)
    pub fn second_derivative_matrix(&self) -> Result<Array2<T>> {
        Self::mat_mul(&self.diff_matrix, &self.diff_matrix)
    }

    /// Apply the second-derivative operator.
    ///
    /// # Errors
    ///
    /// Returns [`cfd_core::error::Error::DimensionMismatch`] when `u` does not
    /// contain exactly one value per collocation point.
    pub fn second_derivative(&self, u: &Array1<T>) -> Result<Array1<T>> {
        let d2 = self.second_derivative_matrix()?;
        Self::mat_vec(&d2, u)
    }

    fn mat_vec(matrix: &Array2<T>, vector: &Array1<T>) -> Result<Array1<T>> {
        let [rows, cols] = matrix.shape();
        if vector.size() != cols {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: cols,
                actual: vector.size(),
            });
        }

        let mut out = Array1::from_elem([rows], scalar::zero::<T>());
        for row in 0..rows {
            let mut sum = scalar::zero::<T>();
            for col in 0..cols {
                sum += matrix[[row, col]] * vector[col];
            }
            out[row] = sum;
        }
        Ok(out)
    }

    fn mat_mul(left: &Array2<T>, right: &Array2<T>) -> Result<Array2<T>> {
        let [left_rows, left_cols] = left.shape();
        let [right_rows, right_cols] = right.shape();
        if left_cols != right_rows {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: left_cols,
                actual: right_rows,
            });
        }

        let mut out = Array2::from_elem([left_rows, right_cols], scalar::zero::<T>());
        for row in 0..left_rows {
            for col in 0..right_cols {
                let mut sum = scalar::zero::<T>();
                for k in 0..left_cols {
                    sum += left[[row, k]] * right[[k, col]];
                }
                out[[row, col]] = sum;
            }
        }
        Ok(out)
    }

    /// Compute Clenshaw-Curtis quadrature weights
    /// Reference: Trefethen (2000), Chapter 12
    pub fn quadrature_weights(&self) -> Result<Vec<T>> {
        let n = self.n;
        if n < 2 {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Need at least 2 points for quadrature".into(),
            ));
        }

        let big_n = n - 1; // N in Trefethen's notation
        let big_n_f64 = big_n as f64;
        let mut w = vec![scalar::zero::<T>(); n];

        let one = scalar::one::<T>();

        // Helper to convert f64 to T
        let to_t = |val: f64| -> T { scalar::from_f64::<T>(val) };

        if big_n.is_multiple_of(2) {
            // N is even
            let val = 1.0 / (big_n_f64 * big_n_f64 - 1.0);
            let w_end = to_t(val);
            w[0] = w_end;
            w[big_n] = w_end;

            for (j, w_j) in w.iter_mut().enumerate().take(big_n).skip(1) {
                let theta_j = PI * (j as f64) / big_n_f64;
                let mut sum = scalar::zero::<T>();

                for k in 1..(big_n / 2) {
                    let coef = 2.0 / (4.0 * (k as f64).powi(2) - 1.0);
                    let term = to_t(coef) * to_t((2.0 * (k as f64) * theta_j).cos());
                    sum += term;
                }

                // Last term for even N
                let last_coef = 1.0 / (big_n_f64 * big_n_f64 - 1.0);
                let last_term = to_t(last_coef) * to_t((big_n_f64 * theta_j).cos());
                sum += last_term;

                *w_j = to_t(2.0 / big_n_f64) * (one - sum);
            }
        } else {
            // N is odd
            let val = 1.0 / (big_n_f64 * big_n_f64);
            let w_end = to_t(val);
            w[0] = w_end;
            w[big_n] = w_end;

            for (j, w_j) in w.iter_mut().enumerate().take(big_n).skip(1) {
                let theta_j = PI * (j as f64) / big_n_f64;
                let mut sum = scalar::zero::<T>();

                for k in 1..=(big_n - 1) / 2 {
                    let coef = 2.0 / (4.0 * (k as f64).powi(2) - 1.0);
                    let term = to_t(coef) * to_t((2.0 * (k as f64) * theta_j).cos());
                    sum += term;
                }

                *w_j = to_t(2.0 / big_n_f64) * (one - sum);
            }
        }

        Ok(w)
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

        let two = scalar::from_f64::<T>(2.0);

        // Barycentric weights for Chebyshev points
        let mut weights = vec![scalar::one::<T>(); self.n];
        weights[0] = scalar::one::<T>() / two;
        weights[self.n - 1] = scalar::from_f64::<T>(if (self.n - 1).is_multiple_of(2) {
            1.0
        } else {
            -1.0
        }) / two;

        for (j, weight) in weights.iter_mut().enumerate().take(self.n - 1).skip(1) {
            *weight = scalar::from_f64::<T>(if j % 2 == 0 { 1.0 } else { -1.0 });
        }

        // Check if x matches a grid point
        for (j, &x_j) in self.points.iter().enumerate() {
            if scalar::abs::<T>(x - x_j) < scalar::from_f64::<T>(1e-14) {
                return Ok(values[j]);
            }
        }

        // Barycentric interpolation
        let mut numer = scalar::zero::<T>();
        let mut denom = scalar::zero::<T>();

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
