//! GMRES solver implementation with Arnoldi iteration and Givens rotations

use super::super::config::IterativeSolverConfig;
use super::super::traits::{Configurable, IterativeLinearSolver, Preconditioner};
use super::{arnoldi, givens};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DMatrix, DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// GMRES(m) solver with restart capability
///
/// Solves non-symmetric linear systems using the Generalized Minimal Residual method.
/// The 'm' parameter controls the maximum Krylov subspace dimension before restart.
///
/// # Algorithm
///
/// 1. Arnoldi process: Build orthonormal basis V for Krylov subspace K_m(A, r0)
/// 2. Modified Gram-Schmidt: Ensure numerical stability of orthogonalization
/// 3. Givens rotations: Solve least-squares problem incrementally
/// 4. Restart: If not converged after m iterations, restart with updated solution
///
/// # Type Parameters
///
/// * `T` - Scalar type (f32 or f64) with real field operations
pub struct GMRES<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
    /// Maximum Krylov subspace dimension before restart
    restart_dim: usize,
}

impl<T: RealField + Copy + FromPrimitive + Debug> GMRES<T> {
    /// Create new GMRES(m) solver
    ///
    /// # Arguments
    ///
    /// * `config` - Solver configuration (tolerance, max iterations)
    /// * `restart_dim` - Maximum Krylov subspace dimension (typically 20-50)
    ///
    /// # Panics
    ///
    /// Panics if restart_dim is zero
    pub fn new(config: IterativeSolverConfig<T>, restart_dim: usize) -> Self {
        assert!(restart_dim > 0, "GMRES restart dimension must be positive");
        Self { config, restart_dim }
    }

    /// Create with default configuration
    ///
    /// Uses restart dimension of 30 (standard for CFD applications)
    #[must_use]
    pub fn default() -> Self {
        Self::new(IterativeSolverConfig::default(), 30)
    }

    /// Solve with preconditioning using GMRES(m) algorithm
    ///
    /// # Arguments
    ///
    /// * `a` - System matrix (n × n sparse CSR)
    /// * `b` - Right-hand side vector (n)
    /// * `preconditioner` - Preconditioner operator
    /// * `x` - Initial guess on entry, solution on exit (n)
    ///
    /// # Errors
    ///
    /// Returns error if:
    /// - Matrix dimensions are incompatible
    /// - Maximum iterations exceeded without convergence
    /// - Numerical breakdown occurs (rare with MGS)
    pub fn solve_preconditioned<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        preconditioner: &P,
        x: &mut DVector<T>,
    ) -> Result<()> {
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

        let m = self.restart_dim.min(n);
        let b_norm = b.norm();
        
        if b_norm < T::from_f64(1e-16).unwrap_or(T::zero()) {
            // RHS is zero, solution is zero
            x.fill(T::zero());
            return Ok(());
        }

        // Workspace allocations (reused across restarts)
        let mut r = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut work = DVector::zeros(n);
        let mut v = DMatrix::zeros(n, m + 1);
        let mut h = DMatrix::zeros(m + 1, m);
        let mut g = DVector::zeros(m + 1);
        let mut y = DVector::zeros(m);
        let mut cs = vec![T::zero(); m];
        let mut sn = vec![T::zero(); m];

        let mut total_iterations = 0;

        // GMRES(m) restart loop
        while total_iterations < self.config.max_iterations {
            // Compute initial residual: r = b - A*x
            arnoldi::spmv(a, x, &mut work);
            r.copy_from(b);
            r -= &work;

            // Apply preconditioner: z = M^{-1} * r
            preconditioner.apply_to(&r, &mut z)?;

            let residual_norm = z.norm();
            
            // Check convergence before starting Arnoldi
            if self.is_converged(residual_norm) {
                tracing::debug!(
                    "GMRES converged in {} total iterations",
                    total_iterations
                );
                return Ok(());
            }

            // Initialize first basis vector: v_0 = z / ||z||
            let inv_norm = T::one() / residual_norm;
            for i in 0..n {
                v[(i, 0)] = z[i] * inv_norm;
            }

            // Initialize RHS for least-squares problem
            g.fill(T::zero());
            g[0] = residual_norm;

            // Arnoldi iteration with Givens rotations
            let mut k = 0;
            while k < m && total_iterations < self.config.max_iterations {
                // Arnoldi step: generate next Krylov basis vector
                match arnoldi::arnoldi_iteration(a, &mut v, &mut h, k, &mut work) {
                    Ok(_norm) => {
                        // Apply Givens rotations to transform H to upper triangular
                        givens::apply_givens_rotation(&mut h, &mut cs, &mut sn, &mut g, k);

                        // Current residual estimate: |g[k+1]|
                        let current_residual = g[k + 1].abs();

                        if self.is_converged(current_residual) {
                            // Converged within this cycle
                            givens::back_substitution(&h, &g, &mut y, k + 1);

                            // Update solution: x = x + V * y
                            for i in 0..n {
                                let mut correction = T::zero();
                                for j in 0..=k {
                                    correction += v[(i, j)] * y[j];
                                }
                                x[i] += correction;
                            }

                            tracing::debug!(
                                "GMRES converged in {} total iterations ({} in last restart)",
                                total_iterations + k + 1,
                                k + 1
                            );
                            return Ok(());
                        }

                        k += 1;
                        total_iterations += 1;
                    }
                    Err(e) => {
                        // Arnoldi breakdown (happy breakdown): residual in Krylov subspace
                        // Compute solution and return
                        if k > 0 {
                            givens::back_substitution(&h, &g, &mut y, k);

                            for i in 0..n {
                                let mut correction = T::zero();
                                for j in 0..k {
                                    correction += v[(i, j)] * y[j];
                                }
                                x[i] += correction;
                            }

                            tracing::debug!(
                                "GMRES converged via breakdown in {} iterations",
                                total_iterations
                            );
                            return Ok(());
                        }
                        return Err(e);
                    }
                }
            }

            // End of restart cycle: update solution with partial result
            givens::back_substitution(&h, &g, &mut y, k);

            for i in 0..n {
                let mut correction = T::zero();
                for j in 0..k {
                    correction += v[(i, j)] * y[j];
                }
                x[i] += correction;
            }

            // Continue to next restart
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Check if residual satisfies convergence criterion
    #[inline]
    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Copy + FromPrimitive + Debug> Configurable<T> for GMRES<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Copy + FromPrimitive + Debug> IterativeLinearSolver<T> for GMRES<T> {
    fn solve<P: Preconditioner<T>>(
        &self,
        a: &CsrMatrix<T>,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<()> {
        match preconditioner {
            Some(p) => self.solve_preconditioned(a, b, p, x),
            None => {
                // Use identity preconditioner
                use super::super::preconditioners::IdentityPreconditioner;
                let identity = IdentityPreconditioner;
                self.solve_preconditioned(a, b, &identity, x)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::super::preconditioners::IdentityPreconditioner;
    use nalgebra_sparse::CooMatrix;

    #[test]
    fn test_gmres_diagonal_matrix() {
        // Simple diagonal matrix: A = diag(1, 2, 3, 4)
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner::default();

        solver.solve_preconditioned(&a, &b, &precond, &mut x).unwrap();

        // Expected solution: x = [1, 2, 3, 4]
        let expected = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let error = (&x - &expected).norm();
        assert!(error < 1e-9, "Solution error: {}", error);
    }

    #[test]
    fn test_gmres_non_symmetric() {
        // Non-symmetric 3x3 matrix
        let n = 3;
        let mut coo = CooMatrix::new(n, n);
        coo.push(0, 0, 4.0);
        coo.push(0, 1, 1.0);
        coo.push(1, 0, 2.0);
        coo.push(1, 1, 5.0);
        coo.push(1, 2, 1.0);
        coo.push(2, 1, 3.0);
        coo.push(2, 2, 6.0);
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![5.0, 8.0, 9.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner::default();

        solver.solve_preconditioned(&a, &b, &precond, &mut x).unwrap();

        // Verify Ax = b
        let ax = &a * &x;
        let residual = (&ax - &b).norm();
        assert!(residual < 1e-9, "Residual: {}", residual);
    }

    #[test]
    fn test_gmres_zero_rhs() {
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::zeros(n);
        let mut x = DVector::from_element(n, 1.0);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner::default();

        solver.solve_preconditioned(&a, &b, &precond, &mut x).unwrap();

        assert!(x.norm() < 1e-14, "Solution should be zero for zero RHS");
    }
}
