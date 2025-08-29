//! BiCGSTAB solver implementation

use super::preconditioners::IdentityPreconditioner;
use super::traits::{LinearSolver, Preconditioner};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
use cfd_core::solvers::LinearSolverConfig;
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// BiCGSTAB solver with efficient memory management
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

    /// Efficient sparse matrix-vector multiplication into pre-allocated buffer
    /// This avoids allocation by reusing the output buffer
    #[inline]
    fn spmv(a: &CsrMatrix<T>, x: &DVector<T>, y: &mut DVector<T>) {
        // Zero out the output vector
        y.fill(T::zero());

        // Perform sparse matrix-vector multiplication manually
        // This is the standard CSR SpMV algorithm
        for i in 0..a.nrows() {
            let row_start = a.row_offsets()[i];
            let row_end = a.row_offsets()[i + 1];

            let mut sum = T::zero();
            for j in row_start..row_end {
                let col_idx = a.col_indices()[j];
                let val = a.values()[j];
                sum += val * x[col_idx];
            }
            y[i] = sum;
        }
    }

    /// Solve with preconditioning using efficient memory management
    ///
    /// # Arguments
    /// * `a` - The coefficient matrix (must be square)
    /// * `b` - The right-hand side vector
    /// * `preconditioner` - The preconditioner to use
    /// * `x` - On entry: initial guess; On exit: solution vector
    ///
    /// # Returns
    /// * `Ok(())` if converged successfully
    /// * `Err(...)` if failed to converge or numerical breakdown occurred
    pub fn solve_preconditioned<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        preconditioner: &P,
        x: &mut DVector<T>, // Changed: mutable reference instead of Option
    ) -> Result<()> {
        // Changed: returns Result<()> instead of Result<DVector<T>>
        let n = b.len();
        if a.nrows() != n || a.ncols() != n {
            return Err(Error::InvalidConfiguration(
                "Matrix dimensions don't match RHS vector".to_string(),
            ));
        }

        if x.len() != n {
            return Err(Error::InvalidConfiguration(
                "Solution vector dimension doesn't match system size".to_string(),
            ));
        }

        // Pre-allocate ALL workspace vectors outside the loop
        let mut r = DVector::zeros(n);
        let mut r0_hat = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n); // Pre-allocated for matrix-vector product
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n); // Pre-allocated for matrix-vector product
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);
        let mut ax = DVector::zeros(n); // Pre-allocated for initial residual

        // Compute initial residual: r = b - A*x
        Self::spmv(a, x, &mut ax); // Allocation-free multiplication
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();

        // Check if already converged
        if self.is_converged(initial_residual_norm) {
            tracing::debug!("BiCGSTAB converged at initial guess");
            return Ok(());
        }

        // Use a robust breakdown tolerance based on machine epsilon
        let epsilon = T::default_epsilon();
        let breakdown_tolerance = epsilon * epsilon; // Simple, robust choice

        r0_hat.copy_from(&r); // Shadow residual

        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for iter in 0..self.config.max_iterations {
            let rho_new = r0_hat.dot(&r);

            // Primary breakdown condition: rho approaching zero
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
            Self::spmv(a, &z, &mut v); // Allocation-free multiplication

            alpha = rho_new / r0_hat.dot(&v);

            // s = r - alpha * v
            s.copy_from(&r);
            s.axpy(-alpha, &v, T::one());

            // REMOVED: Redundant convergence check on s
            // The standard BiCGSTAB algorithm only checks at the end of iteration

            // Solve M*z2 = s and compute t = A*z2
            preconditioner.apply_to(&s, &mut z2)?;
            Self::spmv(a, &z2, &mut t); // Allocation-free multiplication

            let t_dot_t = t.dot(&t);
            if t_dot_t.abs() < breakdown_tolerance {
                // If t is nearly zero, we can't compute omega
                // But we can still update x with what we have
                x.axpy(alpha, &z, T::one());
                // Check if this partial update achieved convergence
                Self::spmv(a, x, &mut ax);
                r.copy_from(b);
                r -= &ax;
                if self.is_converged(r.norm()) {
                    tracing::debug!(
                        "BiCGSTAB converged in {} iterations (t breakdown)",
                        iter + 1
                    );
                    return Ok(());
                }
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            omega = s.dot(&t) / t_dot_t;

            // Update solution: x = x + alpha*z + omega*z2
            x.axpy(alpha, &z, T::one());
            x.axpy(omega, &z2, T::one());

            // Update residual efficiently using swap
            // Instead of r = s - omega * t, we swap and update
            std::mem::swap(&mut r, &mut s); // O(1) swap instead of O(n) copy
            r.axpy(-omega, &t, T::one()); // Now r contains the new residual

            // Single convergence check per iteration (standard BiCGSTAB)
            if self.is_converged(r.norm()) {
                tracing::debug!("BiCGSTAB converged in {} iterations", iter + 1);
                return Ok(());
            }

            // REMOVED: Check on omega - not a standard breakdown condition
            // The primary breakdown is caught by the rho check at loop start

            rho = rho_new;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Solve without preconditioning
    pub fn solve_unpreconditioned(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<()> {
        let identity = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &identity, x)
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
        // For backward compatibility with the trait, we need to allocate here
        // Users should prefer the more efficient solve_preconditioned API
        let mut x = x0.map_or_else(|| DVector::zeros(b.len()), DVector::clone);
        self.solve_unpreconditioned(a, b, &mut x)?;
        Ok(x)
    }

    fn config(&self) -> &LinearSolverConfig<T> {
        &self.config
    }
}
