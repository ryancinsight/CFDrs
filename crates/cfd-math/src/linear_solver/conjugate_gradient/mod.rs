//! Preconditioned Conjugate Gradient (PCG) solver implementation.
//!
//! ## Literature References
//!
//! - Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients for solving linear equations.
//! - Saad, Y. (2003): *Iterative Methods for Sparse Linear Systems*, SIAM. Section 6.7.
//! - Barrett et al. (1994): *Templates for the Solution of Linear Systems*, SIAM.
//! - Golub, G. H., & Van Loan, C. F. (2013): *Matrix Computations* (4th ed.).
//! - Meurant, G. (1999): *Computer Solution of Large Linear Systems*.

use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// Preconditioned Conjugate Gradient solver
pub struct ConjugateGradient<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
}

impl<T: RealField + Copy> ConjugateGradient<T> {
    /// Create new CG solver
    pub const fn new(config: IterativeSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FromPrimitive,
    {
        Self::new(IterativeSolverConfig::default())
    }

    /// Solve with left preconditioning
    pub fn solve_preconditioned<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &DVector<T>,
        preconditioner: &P,
        x: &mut DVector<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        let n = b.len();
        let a_size = a.size();
        if a_size != 0 && a_size != n {
            return Err(Error::InvalidConfiguration(format!(
                "Operator size ({a_size}) doesn't match RHS vector ({n})"
            )));
        }

        let mut r = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut ap = DVector::zeros(n);
        let mut ax = DVector::zeros(n);

        a.apply(x, &mut ax)?;
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        if initial_residual_norm < self.config.tolerance {
            return Ok(monitor);
        }

        preconditioner.apply_to(&r, &mut z)?;
        p.copy_from(&z);

        let epsilon = T::default_epsilon();
        let breakdown_tolerance = epsilon * epsilon;

        let mut r_dot_z = r.dot(&z);
        let r_dot_z_scale = r.norm() * z.norm();
        if r_dot_z.abs() < breakdown_tolerance * (T::one() + r_dot_z_scale) {
            return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
        }
        if r_dot_z < T::zero() {
            return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
        }

        for _iter in 0..self.config.max_iterations {
            a.apply(&p, &mut ap)?;

            let p_dot_ap = p.dot(&ap);
            let p_dot_ap_scale = p.norm() * ap.norm();
            if p_dot_ap.abs() < breakdown_tolerance * (T::one() + p_dot_ap_scale) {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            if p_dot_ap < T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            let alpha = r_dot_z / p_dot_ap;

            x.axpy(alpha, &p, T::one());
            r.axpy(-alpha, &ap, T::one());

            let residual_norm = r.norm();
            monitor.record_residual(residual_norm);

            if residual_norm < self.config.tolerance {
                return Ok(monitor);
            }

            preconditioner.apply_to(&r, &mut z)?;

            let r_dot_z_new = r.dot(&z);
            let r_dot_z_new_scale = r.norm() * z.norm();
            if r_dot_z_new.abs() < breakdown_tolerance * (T::one() + r_dot_z_new_scale) {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            if r_dot_z_new < T::zero() {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            let beta = r_dot_z_new / r_dot_z;

            p *= beta;
            p += &z;

            r_dot_z = r_dot_z_new;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
    }

    /// Solve without preconditioning
    pub fn solve_unpreconditioned<Op: LinearOperator<T> + ?Sized>(
        &self,
        a: &Op,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use super::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }
}

impl<T: RealField + Debug + Copy> Configurable<T> for ConjugateGradient<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy> IterativeLinearSolver<T> for ConjugateGradient<T> {
    fn solve<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &DVector<T>,
        x: &mut DVector<T>,
        preconditioner: Option<&P>,
    ) -> Result<ConvergenceMonitor<T>> {
        if let Some(p) = preconditioner {
            self.solve_preconditioned(a, b, p, x)
        } else {
            self.solve_unpreconditioned(a, b, x)
        }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> super::traits::LinearSolver<T>
    for ConjugateGradient<T>
{
    fn solve_system(
        &self,
        a: &dyn LinearOperator<T>,
        b: &DVector<T>,
        x0: Option<&DVector<T>>,
    ) -> Result<DVector<T>> {
        let mut x = if let Some(initial) = x0 {
            initial.clone()
        } else {
            DVector::zeros(b.len())
        };

        self.solve(
            a,
            b,
            &mut x,
            None::<&super::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::preconditioners::IdentityPreconditioner;
    use super::super::traits::{Configurable, LinearSolver};
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra_sparse::CsrMatrix;

    fn create_simple_spd_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0];

        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let config = IterativeSolverConfig::default();
        let _solver = ConjugateGradient::<f64>::new(config);
    }

    #[test]
    fn test_default_solver() {
        let _solver = ConjugateGradient::<f64>::default();
    }

    #[test]
    fn test_solve_simple_system() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        let ax = &a * &x;
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::from_vec(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_identity_matrix() {
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![1.0, 1.0, 1.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![5.0, 10.0, 15.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        for i in 0..3 {
            assert_relative_eq!(x[i], b[i], epsilon = 1e-10);
        }
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];
        let a = CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![6.0, 9.0, 12.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        assert_relative_eq!(x[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mismatched_dimensions() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0]);
        let mut x = DVector::zeros(2);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::<f64> {
            max_iterations: 1,
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_solve_larger_system() {
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4];
        let values = vec![
            4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0,
        ];
        let a = CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(5);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        let ax = &a * &x;
        for i in 0..5 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_configurable_trait() {
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-8,
            max_iterations: 500,
            ..Default::default()
        };

        let solver = ConjugateGradient::new(config);

        let retrieved_config = solver.config();
        assert_relative_eq!(retrieved_config.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved_config.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_simple_spd_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);

        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());

        let x = result.unwrap();
        let ax = &a * &x;
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }
}
