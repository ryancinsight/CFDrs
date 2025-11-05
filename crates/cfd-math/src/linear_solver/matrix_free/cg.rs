//! Matrix-free Conjugate Gradient solver.
//!
//! This module implements the Conjugate Gradient (CG) method for solving
//! linear systems Ax = b without explicitly storing the matrix A. Instead,
//! A is represented through a LinearOperator that can compute matrix-vector
//! products on demand.
//!
//! CG is optimal for symmetric positive definite systems and converges
//! in at most n iterations for an n x n system.

use super::operator::LinearOperator;
use super::traits::MatrixFreeSolver;
use crate::error::Result;
use crate::linear_solver::config::IterativeSolverConfig;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use cfd_core::error::{Error, ConvergenceErrorKind};

/// Matrix-free Conjugate Gradient solver.
///
/// This solver implements the Conjugate Gradient method for solving
/// symmetric positive definite linear systems using operator abstraction.
/// It requires no matrix storage beyond what the operator needs internally.
pub struct MatrixFreeCG<T: RealField + Copy> {
    /// Solver configuration
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy> MatrixFreeCG<T> {
    /// Create a new matrix-free CG solver with default configuration.
    pub fn new(config: IterativeSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration optimized for symmetric systems.
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        let config = IterativeSolverConfig::default();
        Self::new(config)
    }
}

// Configurable implementation removed for simplicity - will be added back with proper trait bounds

impl<T: RealField + Copy + FromPrimitive> MatrixFreeSolver<T> for MatrixFreeCG<T> {
    /// Solve Ax = b using matrix-free Conjugate Gradient.
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

        // Workspace vectors
        let mut r = vec![T::zero(); n]; // Residual
        let mut p = vec![T::zero(); n]; // Search direction
        let mut ap = vec![T::zero(); n]; // A*p
        let mut z = vec![T::zero(); n]; // Preconditioned residual

        // Initial residual: r = b - A*x
        operator.apply(x, &mut r)?;
        for i in 0..n {
            r[i] = b[i] - r[i];
        }

        // Identity preconditioning for now
        z.copy_from_slice(&r);

        // Initial search direction: p = z
        p.copy_from_slice(&z);

        // Initial residual norms
        let mut _r_norm_sq = self.dot(&r, &r);
        let mut z_norm_sq = self.dot(&z, &z);

        let mut iterations = 0;
        let mut converged = false;

        while iterations < self.config.max_iterations && !converged {
            // Compute A*p
            operator.apply(&p, &mut ap)?;

            // Compute step size: alpha = (r*z) / (p*A*p)
            let p_ap = self.dot(&p, &ap);
            if p_ap == T::zero() {
                break; // Breakdown
            }
            let alpha = z_norm_sq / p_ap;

            // Update solution: x = x + alpha*p
            for i in 0..n {
                x[i] += alpha * p[i];
            }

            // Update residual: r = r - alpha*A*p
            for i in 0..n {
                r[i] -= alpha * ap[i];
            }

            // Identity preconditioning for now
            z.copy_from_slice(&r);

            // Check convergence
            let new_r_norm_sq = self.dot(&r, &r);
            if new_r_norm_sq.sqrt() < self.config.tolerance {
                converged = true;
                break;
            }

            // Compute beta = (new_r*z) / (old_r*z)
            let new_z_norm_sq = self.dot(&z, &z);
            let beta = new_z_norm_sq / z_norm_sq;

            // Update search direction: p = z + beta*p
            for i in 0..n {
                p[i] = z[i] + beta * p[i];
            }

            // Update norms for next iteration
            _r_norm_sq = new_r_norm_sq;
            z_norm_sq = new_z_norm_sq;

            iterations += 1;
        }

        if !converged && iterations >= self.config.max_iterations {
            return Err(Error::Convergence(
                ConvergenceErrorKind::MaxIterationsExceeded { max: self.config.max_iterations },
            ));
        }

        Ok(())
    }
}

impl<T: RealField + Copy> MatrixFreeCG<T> {
    /// Compute dot product of two vectors.
    fn dot(&self, a: &[T], b: &[T]) -> T {
        debug_assert_eq!(a.len(), b.len());
        let mut sum = T::zero();
        for i in 0..a.len() {
            sum += a[i] * b[i];
        }
        sum
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::operator::IdentityOperator;
    use approx::assert_relative_eq;

    #[test]
    fn test_cg_identity_operator() {
        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = MatrixFreeCG::new(config);
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
    fn test_cg_convergence() {
        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(100);
        let solver = MatrixFreeCG::new(config);

        // Simple 2x2 SPD system: [[4, 1], [1, 3]] [x, y] = [1, 2]
        // Solution: x = 1/11, y = 7/11
        struct SimpleOperator;
        impl LinearOperator<f64> for SimpleOperator {
            fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
                y[0] = 4.0 * x[0] + 1.0 * x[1];
                y[1] = 1.0 * x[0] + 3.0 * x[1];
                Ok(())
            }

            fn size(&self) -> usize { 2 }

            fn is_symmetric(&self) -> bool { true }

            fn is_positive_definite(&self) -> Option<bool> { Some(true) }
        }

        let operator = SimpleOperator;
        let b = vec![1.0, 2.0];
        let mut x = vec![0.0; 2];

        solver.solve(&operator, &b, &mut x).unwrap();

        assert_relative_eq!(x[0], 1.0 / 11.0, epsilon = 1e-8);
        assert_relative_eq!(x[1], 7.0 / 11.0, epsilon = 1e-8);
    }
}
