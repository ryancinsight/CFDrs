//! Preconditioned Conjugate Gradient solver implementation

use super::traits::{LinearSolver, Preconditioner};
use crate::vector_ops::SimdVectorOps;
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use cfd_core::solvers::LinearSolverConfig;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// Preconditioned Conjugate Gradient solver with memory management
pub struct ConjugateGradient<T: RealField + Copy> {
    config: LinearSolverConfig<T>,
}

impl<T: RealField + Copy> ConjugateGradient<T> {
    /// Create new CG solver
    pub const fn new(config: LinearSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        Self::new(LinearSolverConfig::default())
    }

    /// Solve with preconditioning and memory management
    pub fn solve_preconditioned<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        preconditioner: &P,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        // Pre-allocate all workspace vectors to avoid allocations in the loop
        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);
        let mut r = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut ap;

        // Compute initial residual: r = b - A*x
        let ax = a * &x;
        r.copy_from(b);
        r -= &ax;

        // Apply preconditioner: M*z = r
        preconditioner.apply_to(&r, &mut z)?;
        p.copy_from(&z);

        let mut rzold = r.dot(&z);

        // PCG iterations with in-place operations
        for iter in 0..self.config.max_iterations {
            // Compute A*p
            ap = a * &p;

            let alpha = rzold / p.dot(&ap);

            // Update solution: x = x + alpha * p
            x.axpy(alpha, &p, T::one());

            // Update residual: r = r - alpha * ap
            r.axpy(-alpha, &ap, T::one());

            let residual_norm = r.norm();
            if self.is_converged(residual_norm) {
                tracing::debug!("PCG converged in {} iterations", iter + 1);
                return Ok(x);
            }

            // Apply preconditioner: M*z = r
            preconditioner.apply_to(&r, &mut z)?;

            let rznew = r.dot(&z);
            let beta = rznew / rzold;

            // Update search direction: p = z + beta * p
            p *= beta;
            p += &z;

            rzold = rznew;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy + FromPrimitive + Send + Sync> LinearSolver<T>
    for ConjugateGradient<T>
{
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidInput(
                "Matrix dimensions must match vector size".to_string(),
            ));
        }

        let mut x = x0.map(|v| v.clone()).unwrap_or_else(|| DVector::zeros(n));
        let mut r = b - a * &x;

        // Use SIMD operations for large vectors
        let mut r_norm_sq = if n > 1000 {
            r.simd_norm().powi(2)
        } else {
            r.norm_squared()
        };

        if r_norm_sq < self.config.tolerance * self.config.tolerance {
            return Ok(x);
        }

        let mut p = r.clone();

        for _ in 0..self.config.max_iterations {
            let ap = a * &p;

            // Use SIMD dot product for large vectors
            let pap = if n > 1000 {
                p.simd_dot(&ap)
            } else {
                p.dot(&ap)
            };

            if pap.abs() < T::from_f64(1e-14).unwrap_or_else(|| T::zero()) {
                break;
            }

            let alpha = r_norm_sq / pap;

            // Use axpy for in-place updates
            for i in 0..n {
                x[i] = x[i] + alpha * p[i];
                r[i] = r[i] - alpha * ap[i];
            }

            // Use SIMD norm for convergence check
            let r_norm_sq_new = if n > 1000 {
                r.simd_norm().powi(2)
            } else {
                r.norm_squared()
            };

            if r_norm_sq_new < self.config.tolerance * self.config.tolerance {
                return Ok(x);
            }

            let beta = r_norm_sq_new / r_norm_sq;

            // Update p = r + beta * p
            for i in 0..n {
                p[i] = r[i] + beta * p[i];
            }

            r_norm_sq = r_norm_sq_new;
        }

        Ok(x)
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}
