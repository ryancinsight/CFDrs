//! BiCGSTAB solver implementation

use super::traits::{LinearSolver, Preconditioner};
use super::preconditioners::IdentityPreconditioner;
use cfd_core::error::{Error, Result, NumericalErrorKind, ConvergenceErrorKind};
use cfd_core::solver::LinearSolverConfig;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// BiCGSTAB solver with optimized memory management
pub struct BiCGSTAB<T: RealField + Copy> {
    config: LinearSolverConfig<T>,
}

impl<T: RealField + Copy> BiCGSTAB<T> {
    /// Create new BiCGSTAB solver
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

    /// Solve with preconditioning and optimized memory management
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

        // Pre-allocate workspace vectors
        let mut x = x0.map_or_else(|| DVector::zeros(n), DVector::clone);
        let mut r = DVector::zeros(n);
        let mut r0_hat = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n);
        let mut s = DVector::zeros(n);
        let mut t;
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);

        // Compute initial residual: r = b - A*x
        let ax = a * &x;
        r.copy_from(b);
        r -= &ax;
        
        let initial_residual_norm = r.norm();
        // Use a more robust breakdown tolerance
        // The tolerance should be at least sqrt(epsilon) to avoid false positives
        let epsilon = T::default_epsilon();
        let sqrt_epsilon = epsilon.sqrt();
        let breakdown_tolerance = (initial_residual_norm * sqrt_epsilon).max(epsilon);
        
        r0_hat.copy_from(&r); // Shadow residual
        
        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for iter in 0..self.config.max_iterations {
            let rho_new = r0_hat.dot(&r);
            
            if rho_new.abs() < breakdown_tolerance {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            let beta = (rho_new / rho) * (alpha / omega);
            
            // p = r + beta * (p - omega * v) - using in-place operations
            p.axpy(-omega, &v, T::one());
            p *= beta;
            p += &r;

            // Solve M*z = p and compute v = A*z
            preconditioner.apply_to(&p, &mut z)?;
            v = a * &z;
            
            alpha = rho_new / r0_hat.dot(&v);
            
            // s = r - alpha * v
            s.copy_from(&r);
            s.axpy(-alpha, &v, T::one());
            
            if self.is_converged(s.norm()) {
                x.axpy(alpha, &z, T::one());
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }

            // Solve M*z2 = s and compute t = A*z2
            preconditioner.apply_to(&s, &mut z2)?;
            t = a * &z2;
            
            omega = s.dot(&t) / t.dot(&t);
            
            // Update solution: x = x + alpha*z + omega*z2
            x.axpy(alpha, &z, T::one());
            x.axpy(omega, &z2, T::one());
            
            // Update residual: r = s - omega * t
            r.copy_from(&s);
            r.axpy(-omega, &t, T::one());
            
            if self.is_converged(r.norm()) {
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(x);
            }
            
            if omega.abs() < breakdown_tolerance {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            rho = rho_new;
        }

        Err(Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { 
            max: self.config.max_iterations 
        }))
    }

    /// Check if residual satisfies convergence criteria
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy> LinearSolver<T> for BiCGSTAB<T> {
    fn solve(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        // Use identity preconditioner for unpreconditioned solve
        let identity = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &identity, x0)
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}