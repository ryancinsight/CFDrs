//! Matrix-free GMRES solver.
//!
//! This module implements the Generalized Minimal Residual (GMRES) method
//! for solving linear systems Ax = b using operator abstraction. GMRES is
//! suitable for general (non-symmetric) linear systems and can handle
//! indefinite and singular matrices.
//!
//! The implementation uses restarted GMRES(m) with Givens rotations for
//! orthogonalization, providing both stability and efficiency.

use super::operator::LinearOperator;
use super::traits::MatrixFreeSolver;
use crate::error::Result;
use crate::linear_solver::config::IterativeSolverConfig;
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::FromPrimitive;
use cfd_core::error::{Error, ConvergenceErrorKind};

/// Matrix-free GMRES solver.
///
/// This solver implements the restarted GMRES(m) method for solving
/// general linear systems using operator abstraction. It maintains
/// a Krylov subspace of dimension m and restarts when this limit
/// is reached.
pub struct MatrixFreeGMRES<T: RealField + Copy> {
    /// Solver configuration
    config: IterativeSolverConfig<T>,
    /// Krylov subspace dimension (restart parameter)
    restart_dim: usize,
}

impl<T: RealField + Copy> MatrixFreeGMRES<T> {
    /// Create a new matrix-free GMRES solver.
    ///
    /// # Arguments
    ///
    /// * `config` - Solver configuration
    /// * `restart_dim` - Krylov subspace dimension (m in GMRES(m))
    pub fn new(config: IterativeSolverConfig<T>, restart_dim: usize) -> Self {
        Self { config, restart_dim }
    }

    /// Create with default configuration and restart dimension.
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        let config = IterativeSolverConfig::default();
        Self::new(config, 30) // Standard restart dimension
    }
}

// Configurable implementation removed for simplicity - will be added back with proper trait bounds

impl<T: RealField + Copy + FromPrimitive> MatrixFreeSolver<T> for MatrixFreeGMRES<T> {
    /// Solve Ax = b using matrix-free restarted GMRES.
    ///
    /// This simplified implementation uses identity preconditioning.
    /// Full preconditioner support will be added in future iterations.
    ///
    /// # Arguments
    ///
    /// * `operator` - Linear operator representing A
    /// * `b` - Right-hand side vector
    /// * `x` - Initial guess (input) / solution (output)
    ///
    /// # Errors
    ///
    /// Returns error if convergence fails or operator application fails.
    fn solve(&self, operator: &dyn LinearOperator<T>, b: &[T], x: &mut [T]) -> Result<()> {
        let n = operator.size();
        if b.len() != n || x.len() != n {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let m = self.restart_dim;
        if m == 0 {
            return Err(Error::InvalidConfiguration(
                "Restart dimension must be positive".to_string(),
            ));
        }

        // GMRES(m) restart loop
        let mut total_iterations = 0;
        let mut converged = false;

        while total_iterations < self.config.max_iterations && !converged {
            // Solve one GMRES(m) cycle
            let (iterations_used, residual_norm) = self.gmres_cycle(operator, b, x, m)?;

            total_iterations += iterations_used;

            // Check convergence
            if residual_norm < self.config.tolerance {
                converged = true;
            }
        }

        if !converged && total_iterations >= self.config.max_iterations {
            return Err(Error::Convergence(
                ConvergenceErrorKind::MaxIterationsExceeded {
                    max: self.config.max_iterations
                },
            ));
        }

        Ok(())
    }
}

impl<T: RealField + Copy + FromPrimitive> MatrixFreeGMRES<T> {
    /// Perform one GMRES(m) cycle.
    fn gmres_cycle(
        &self,
        operator: &dyn LinearOperator<T>,
        b: &[T],
        x: &mut [T],
        m: usize,
    ) -> Result<(usize, T)> {
        let n = operator.size();

        // Arnoldi workspace
        let mut v = DMatrix::zeros(n, m + 1); // Krylov basis vectors
        let mut h = DMatrix::zeros(m + 1, m); // Hessenberg matrix
        let mut g = DVector::zeros(m + 1); // Right-hand side for least squares

        // Givens rotation storage
        let mut cs = vec![T::zero(); m];
        let mut sn = vec![T::zero(); m];

        // Temporary vectors
        let mut w = vec![T::zero(); n];
        let mut r = vec![T::zero(); n];

        // Initial residual: r = b - A*x
        operator.apply(x, &mut r)?;
        for i in 0..n {
            r[i] = b[i] - r[i];
        }

        // Identity preconditioning for now

        // Initial residual norm
        let r_norm = self.vector_norm(&r);
        if r_norm == T::zero() {
            return Ok((0, T::zero())); // Exact solution
        }

        // First basis vector: v1 = r / ||r||
        let r_norm_inv = T::one() / r_norm;
        for i in 0..n {
            v[(i, 0)] = r[i] * r_norm_inv;
        }
        g[0] = r_norm;

        let mut j = 0;
        while j < m {
            // Arnoldi iteration: w = A*v_j
            let vj_col = v.column(j);
            operator.apply(vj_col.as_slice(), &mut w)?;

            // Identity preconditioning for now

            // Modified Gram-Schmidt orthogonalization
            for i in 0..=j {
                let vi_col = v.column(i);
                let h_ij = self.dot_product(vi_col.as_slice(), &w);
                h[(i, j)] = h_ij;

                // w = w - h_ij * v_i
                for k in 0..n {
                    w[k] -= h_ij * vi_col[k];
                }
            }

            // Compute norm of w
            let w_norm = self.vector_norm(&w);
            h[(j + 1, j)] = w_norm;

            if w_norm == T::zero() {
                break; // Breakdown
            }

            // Normalize: v_{j+1} = w / ||w||
            let w_norm_inv = T::one() / w_norm;
            for i in 0..n {
                v[(i, j + 1)] = w[i] * w_norm_inv;
            }

            // Apply Givens rotations to Hessenberg matrix
            for i in 0..j {
                let temp = cs[i] * h[(i, j)] + sn[i] * h[(i + 1, j)];
                h[(i + 1, j)] = -sn[i] * h[(i, j)] + cs[i] * h[(i + 1, j)];
                h[(i, j)] = temp;
            }

            // Compute new Givens rotation
            let (cs_j, sn_j) = self.givens_rotation(h[(j, j)], h[(j + 1, j)]);
            cs[j] = cs_j;
            sn[j] = sn_j;

            // Apply rotation to Hessenberg matrix
            h[(j, j)] = cs_j * h[(j, j)] + sn_j * h[(j + 1, j)];
            h[(j + 1, j)] = T::zero();

            // Apply rotation to right-hand side
            let temp = cs_j * g[j] + sn_j * g[j + 1];
            g[j + 1] = -sn_j * g[j] + cs_j * g[j + 1];
            g[j] = temp;

            j += 1;

            // Check convergence
            let residual = g[j].abs();
            if residual < self.config.tolerance {
                break;
            }
        }

        // Solve upper triangular system: H*y = g
        let mut y = DVector::zeros(j);
        for i in (0..j).rev() {
            let mut sum = T::zero();
            for k in (i + 1)..j {
                sum += h[(i, k)] * y[k];
            }
            y[i] = (g[i] - sum) / h[(i, i)];
        }

        // Update solution: x = x + V*y
        for i in 0..n {
            let mut update = T::zero();
            for k in 0..j {
                update += v[(i, k)] * y[k];
            }
            x[i] += update;
        }

        let final_residual = g[j].abs();
        Ok((j, final_residual))
    }

    /// Compute Givens rotation parameters.
    fn givens_rotation(&self, a: T, b: T) -> (T, T) {
        let r = (a * a + b * b).sqrt();
        if r == T::zero() {
            (T::one(), T::zero())
        } else {
            (a / r, b / r)
        }
    }

    /// Compute dot product of two vectors.
    fn dot_product(&self, a: &[T], b: &[T]) -> T {
        debug_assert_eq!(a.len(), b.len());
        let mut sum = T::zero();
        for i in 0..a.len() {
            sum += a[i] * b[i];
        }
        sum
    }

    /// Compute L2 norm of a vector.
    fn vector_norm(&self, v: &[T]) -> T {
        self.dot_product(v, v).sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::operator::IdentityOperator;
    use approx::assert_relative_eq;

    #[test]
    fn test_gmres_identity_operator() {
        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = MatrixFreeGMRES::new(config, 10);
        let operator = IdentityOperator::new(5);

        let b = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let mut x = vec![0.0; 5];

        solver.solve(&operator, &b, &mut x).unwrap();

        // For identity matrix, solution should equal RHS
        for i in 0..5 {
            assert_relative_eq!(x[i], b[i], epsilon = 1e-8);
        }
    }

    #[test]
    fn test_gmres_simple_system() {
        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(100);
        let solver = MatrixFreeGMRES::new(config, 5);

        // Simple upper triangular system
        struct UpperTriangularOperator;
        impl LinearOperator<f64> for UpperTriangularOperator {
            fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
                y[0] = 2.0 * x[0] + x[1];
                y[1] = 3.0 * x[1];
                Ok(())
            }

            fn size(&self) -> usize { 2 }
        }

        let operator = UpperTriangularOperator;
        let b = vec![5.0, 6.0]; // Solution should be [1, 2]
        let mut x = vec![0.0; 2];

        solver.solve(&operator, &b, &mut x).unwrap();

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-6);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-6);
    }
}
