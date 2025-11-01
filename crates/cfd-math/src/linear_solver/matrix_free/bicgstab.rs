//! Matrix-free BiCGSTAB solver.
//!
//! This module implements the Bi-Conjugate Gradient Stabilized (BiCGSTAB)
//! method for solving linear systems Ax = b using operator abstraction.
//! BiCGSTAB is a robust solver for general non-symmetric linear systems,
//! particularly effective for CFD applications with complex physics.

use super::operator::LinearOperator;
use super::traits::MatrixFreeSolver;
use crate::linear_solver::config::IterativeSolverConfig;
use cfd_core::error::{Error, Result, ConvergenceErrorKind, NumericalErrorKind};
use nalgebra::RealField;

/// Matrix-free BiCGSTAB solver.
///
/// This solver implements the Bi-Conjugate Gradient Stabilized method
/// for solving general (non-symmetric) linear systems using operator
/// abstraction. BiCGSTAB typically converges faster than GMRES for
/// many CFD applications while being more memory-efficient.
pub struct MatrixFreeBiCGSTAB<T> 
where
    T: RealField + Copy + Default + From<f64> + 'static,
{
    /// Solver configuration
    config: IterativeSolverConfig<T>,
}

impl<T> MatrixFreeBiCGSTAB<T> 
where
    T: RealField + Copy + Default + From<f64> + 'static,
{
    /// Create a new matrix-free BiCGSTAB solver.
    pub fn new(config: IterativeSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration.
    pub fn default() -> Self
    where
        T: From<f64> + Default,
    {
        let config = IterativeSolverConfig::default();
        Self::new(config)
    }
}

use crate::linear_solver::traits::Configurable;

impl<T> Configurable<T> for MatrixFreeBiCGSTAB<T> 
where
    T: RealField + Copy + Default + From<f64> + 'static,
{
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T> MatrixFreeSolver<T> for MatrixFreeBiCGSTAB<T> 
where
    T: RealField + Copy + Default + From<f64> + std::fmt::Debug + 'static,
{
    /// Solve Ax = b using matrix-free BiCGSTAB.
    ///
    /// This implementation uses identity preconditioning. Full preconditioner
    /// support will be added in future iterations.
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

        // BiCGSTAB workspace vectors
        let mut r = vec![T::zero(); n];     // Residual vector
        let mut r_hat = vec![T::zero(); n]; // Shadow residual
        let mut p = vec![T::zero(); n];     // Search direction
        let mut p_hat = vec![T::zero(); n]; // Shadow search direction
        let mut v = vec![T::zero(); n];     // A * p
        let mut s = vec![T::zero(); n];     // Intermediate vector
        let mut t = vec![T::zero(); n];     // A * s

        // Initial residual: r = b - A*x
        operator.apply(x, &mut r)?;
        for i in 0..n {
            r[i] = b[i] - r[i];
        }

        let mut r_norm = self.vector_norm(&r);
        if r_norm < self.config.tolerance {
            return Ok(());
        }

        // Initialize shadow residual: r̂ = r (choose arbitrary)
        r_hat.copy_from_slice(&r);

        // Initialize scalars
        let mut rho_old = self.dot_product(&r_hat, &r);
        let mut _alpha = T::from_f64(1.0).unwrap();
        let mut _omega = T::from_f64(1.0).unwrap();

        // Initial search directions
        p.copy_from_slice(&r);
        p_hat.copy_from_slice(&r);

        let mut iteration = 0;

        while iteration < self.config.max_iterations {
            // v = A * p̂
            operator.apply(&p_hat, &mut v)?;

            // alpha = rho_old / (r̂·v)
            let rho_new = self.dot_product(&r_hat, &v);
            if rho_new == T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }
            _alpha = rho_old / rho_new;

            // s = r - alpha * v
            for i in 0..n {
                s[i] = r[i] - _alpha * v[i];
            }

            let s_norm = self.vector_norm(&s);
            if s_norm < self.config.tolerance {
                // Early convergence
                for i in 0..n {
                    x[i] += _alpha * p[i];
                }
                return Ok(());
            }

            // t = A * s
            operator.apply(&s, &mut t)?;

            // omega = (t·s) / (t·t)
            let ts_dot = self.dot_product(&t, &s);
            let tt_dot = self.dot_product(&t, &t);

            if tt_dot == T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            _omega = ts_dot / tt_dot;

            // x = x + alpha * p + omega * s
            for i in 0..n {
                x[i] += _alpha * p[i] + _omega * s[i];
            }

            // r = s - omega * t
            for i in 0..n {
                r[i] = s[i] - _omega * t[i];
            }

            r_norm = self.vector_norm(&r);
            if r_norm < self.config.tolerance {
                return Ok(());
            }

            // rho_old = r̂·r
            rho_old = self.dot_product(&r_hat, &r);

            if rho_old == T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            // beta = (rho_old / rho_new) * (alpha / omega)
            let beta = (rho_old / rho_new) * (_alpha / _omega);

            // p = r + beta * (p - omega * v)
            // p̂ = r̂ + beta * (p̂ - omega * v̂)  (but v̂ is approximated by v)
            for i in 0..n {
                p[i] = r[i] + beta * (p[i] - _omega * v[i]);
                p_hat[i] = r_hat[i] + beta * (p_hat[i] - _omega * v[i]);
            }

            iteration += 1;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }
}

impl<T> MatrixFreeBiCGSTAB<T> 
where
    T: RealField + Copy + Default + From<f64> + 'static,
{
    /// Compute dot product of two vectors.
    fn dot_product(&self, a: &[T], b: &[T]) -> T {
        debug_assert_eq!(a.len(), b.len());
        let mut sum = T::zero();
        for i in 0..a.len() {
            sum = sum + a[i] * b[i];
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
    fn test_bicgstab_identity_operator() {
        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = MatrixFreeBiCGSTAB::new(config);
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
    fn test_bicgstab_simple_system() {
        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(100);
        let solver = MatrixFreeBiCGSTAB::new(config);

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

    #[test]
    fn test_bicgstab_diagonal_system() {
        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = MatrixFreeBiCGSTAB::new(config);

        // Diagonal system: Ax = b where A is diagonal matrix
        struct DiagonalOperator;
        impl LinearOperator<f64> for DiagonalOperator {
            fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<()> {
                y[0] = x[0];
                y[1] = 2.0 * x[1];
                y[2] = x[2] / 2.0;
                Ok(())
            }

            fn size(&self) -> usize { 3 }
        }

        let operator = DiagonalOperator;
        let b = vec![1.0, 4.0, 1.5]; // Solution should be [1, 2, 3]
        let mut x = vec![0.0; 3];

        solver.solve(&operator, &b, &mut x).unwrap();

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-8);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-8);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-8);
    }
}
