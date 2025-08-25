//! Chebyshev polynomial operations for spectral methods
//!
//! Reference: Trefethen, L.N. (2000). "Spectral Methods in MATLAB"

use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;
use cfd_core::error::Result;
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
            let x = T::from_f64(theta.cos())
                .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                    "Cannot convert collocation point".into()
                ))?;
            points.push(x);
        }
        
        Ok(points)
    }
    
    /// Construct Chebyshev differentiation matrix
    /// Based on Trefethen (2000), Chapter 6
    fn differentiation_matrix(n: usize, points: &[T]) -> Result<DMatrix<T>> {
        let mut d = DMatrix::zeros(n, n);
        let two = T::from_f64(2.0)
            .ok_or_else(|| cfd_core::error::Error::InvalidConfiguration(
                "Cannot convert constant".into()
            ))?;
        
        // c_i = 2 for i = 0 or N, 1 otherwise
        let mut c = vec![T::one(); n];
        c[0] = two;
        c[n - 1] = two;
        
        // Fill differentiation matrix
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    let num = c[i] * T::from_i32(if (i + j) % 2 == 0 { 1 } else { -1 }).expect("CRITICAL: Add proper error handling");
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
    #[must_use] pub fn differentiate(&self, u: &DVector<T>) -> DVector<T> {
        &self.diff_matrix * u
    }
    
    /// Get collocation points
    #[must_use] pub fn points(&self) -> &[T] {
        &self.points
    }
    
    /// Get second derivative matrix (D²)
    pub fn second_derivative_matrix(&self) -> Result<DMatrix<T>> {
        // Second derivative is D * D
        Ok(&self.diff_matrix * &self.diff_matrix)
    }
}

/// Chebyshev differentiation operations
pub struct ChebyshevDifferentiation<T: RealField + Copy> {
    basis: ChebyshevPolynomial<T>,
}

impl<T: RealField + FromPrimitive + Copy> ChebyshevDifferentiation<T> {
    pub fn new(n: usize) -> Result<Self> {
        Ok(Self {
            basis: ChebyshevPolynomial::new(n)?,
        })
    }
    
    /// Compute first derivative
    #[must_use] pub fn first_derivative(&self, u: &DVector<T>) -> DVector<T> {
        self.basis.differentiate(u)
    }
    
    /// Compute second derivative
    #[must_use] pub fn second_derivative(&self, u: &DVector<T>) -> DVector<T> {
        let du = self.basis.differentiate(u);
        self.basis.differentiate(&du)
    }
}