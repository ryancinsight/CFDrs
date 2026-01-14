//! BiCGSTAB solver implementation.

use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::fmt::Debug;

/// BiCGSTAB solver with efficient memory management
///
/// # Convergence Theory
///
/// ## Acceleration over BiCG
///
/// BiCGSTAB (BiConjugate Gradient Stabilized) addresses instability issues in the
/// standard BiCG method by introducing additional stabilization terms. The algorithm
/// combines BiCG with GMRES-like stabilization to prevent divergence.
///
/// Unlike BiCG, which may exhibit erratic convergence or breakdown, BiCGSTAB
/// provides smoother convergence behavior and better numerical stability.
///
/// ## Convergence Properties
///
/// BiCGSTAB converges for non-symmetric matrices where the field of values
/// lies in the right half-plane. The convergence rate depends on:
/// - Condition number of the system matrix
/// - Effectiveness of preconditioning
/// - Initial residual smoothness
///
/// ## Breakdown Prevention
///
/// BiCGSTAB handles two types of breakdown:
/// 1. **Primary breakdown**: ρ_new = 0 (division by zero)
/// 2. **Secondary breakdown**: α = 0 or ω = 0 (stagnation)
///
/// Robust implementations include checks for near-breakdown conditions.
///
/// # References
///
/// - Van der Vorst, H. A. (1992). Bi-CGSTAB: A fast and smoothly converging variant
///   of Bi-CG for the solution of nonsymmetric linear systems. *SIAM Journal on
///   Scientific and Statistical Computing*, 13(2), 631-644.
///   See Algorithm 1 and convergence analysis in Section 3.
/// - Sleijpen, G. L. G., & Fokkema, D. R. (1993). BiCGstab(l) and other hybrid
///   methods for systems of nonsymmetric linear equations. *Numerical Algorithms*,
///   7(3), 347-369.
pub struct BiCGSTAB<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
}

fn is_finite_scalar<T: RealField + Copy>(x: T) -> bool {
    x.to_subset().is_none_or(f64::is_finite)
}

impl<T: RealField + Copy> BiCGSTAB<T> {
    /// Create new BiCGSTAB solver
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
        let mut r0_hat = DVector::zeros(n);
        let mut p = DVector::zeros(n);
        let mut v = DVector::zeros(n);
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);
        let mut ax = DVector::zeros(n);

        a.apply(x, &mut ax)?;
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        if self.is_converged(initial_residual_norm) {
            tracing::debug!("BiCGSTAB converged at initial guess");
            return Ok(monitor);
        }

        let epsilon = T::default_epsilon();
        let breakdown_tolerance = epsilon * epsilon;

        r0_hat.copy_from(&r);

        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for _iter in 0..self.config.max_iterations {
            let rho_new = r0_hat.dot(&r);
            let rho_scale = r0_hat.norm() * r.norm();

            if rho_new.abs() < breakdown_tolerance * (T::one() + rho_scale) {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }

            if omega.abs() < breakdown_tolerance * (T::one() + alpha.abs()) {
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: r.norm().to_subset().unwrap_or(0.0),
                    },
                ));
            }

            let beta = (rho_new / rho) * (alpha / omega);

            p.axpy(-omega, &v, T::one());
            p *= beta;
            p += &r;

            preconditioner.apply_to(&p, &mut z)?;
            a.apply(&z, &mut v)?;

            let r0_hat_dot_v = r0_hat.dot(&v);
            let r0_hat_dot_v_scale = r0_hat.norm() * v.norm();
            if r0_hat_dot_v.abs() < breakdown_tolerance * (T::one() + r0_hat_dot_v_scale) {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            alpha = rho_new / r0_hat_dot_v;
            if !is_finite_scalar(alpha) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }

            s.copy_from(&r);
            s.axpy(-alpha, &v, T::one());

            let s_norm = s.norm();
            if self.is_converged(s_norm) {
                x.axpy(alpha, &z, T::one());
                monitor.record_residual(s_norm);
                return Ok(monitor);
            }

            preconditioner.apply_to(&s, &mut z2)?;
            a.apply(&z2, &mut t)?;

            let t_dot_t = t.dot(&t);
            let t_dot_t_scale = t.norm() * t.norm();
            if t_dot_t.abs() < breakdown_tolerance * (T::one() + t_dot_t_scale) {
                x.axpy(alpha, &z, T::one());
                a.apply(x, &mut ax)?;
                r.copy_from(b);
                r -= &ax;
                let final_norm = r.norm();
                monitor.record_residual(final_norm);
                if self.is_converged(final_norm) {
                    return Ok(monitor);
                }
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: final_norm.to_subset().unwrap_or(0.0),
                    },
                ));
            }

            omega = t.dot(&s) / t_dot_t;
            if !is_finite_scalar(omega) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }
            x.axpy(alpha, &z, T::one());
            x.axpy(omega, &z2, T::one());

            r.copy_from(&s);
            r.axpy(-omega, &t, T::one());

            let residual_norm = r.norm();
            if !is_finite_scalar(residual_norm) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }
            monitor.record_residual(residual_norm);

            if self.is_converged(residual_norm) {
                return Ok(monitor);
            }

            if omega.abs() < breakdown_tolerance * (T::one() + residual_norm) {
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: residual_norm.to_subset().unwrap_or(0.0),
                    },
                ));
            }

            rho = rho_new;
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

    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy> Configurable<T> for BiCGSTAB<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy> IterativeLinearSolver<T> for BiCGSTAB<T> {
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
    for BiCGSTAB<T>
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

    fn create_nonsymmetric_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![5.0, 1.0, 2.0, 4.0, 1.0, 2.0, 3.0];

        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    fn create_diagonal_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];

        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    #[test]
    fn test_new_solver() {
        let config = IterativeSolverConfig::default();
        let _solver = BiCGSTAB::<f64>::new(config);
    }

    #[test]
    fn test_default_solver() {
        let _solver = BiCGSTAB::<f64>::default();
    }

    #[test]
    fn test_solve_simple_system() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
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
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::from_vec(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![2.0, 3.0, 4.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_mismatched_dimensions() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![1.0, 2.0]);
        let mut x = DVector::zeros(2);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::<f64> {
            max_iterations: 1,
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_configurable_trait() {
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-8,
            max_iterations: 500,
            ..Default::default()
        };

        let solver = BiCGSTAB::new(config);

        let retrieved_config = solver.config();
        assert_relative_eq!(retrieved_config.tolerance, 1e-8, epsilon = 1e-10);
        assert_eq!(retrieved_config.max_iterations, 500);
    }

    #[test]
    fn test_linear_solver_trait() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);

        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());

        let x = result.unwrap();
        let ax = &a * &x;
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_already_converged() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![2.0, 3.0, 4.0]);
        let mut x = DVector::from_vec(vec![1.0, 1.0, 1.0]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 1.0, epsilon = 1e-10);
    }
}
