//! BiCGSTAB solver implementation.

use super::array_ops::{
    assign_residual, axpy, copy_array, dot, norm, scale_add, validate_vector_len, vector_len,
};
use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
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
/// # Theorem (BiCGSTAB Residual Reduction)
///
/// Each BiCGSTAB iteration applies a degree-2 polynomial to the residual:
/// $r_k = p_{2k}(A)\,r_0$ where $p_{2k} = \prod_{j=1}^{k}(I - \omega_j A)(I - \alpha_j A) \cdot p_{2(k-1)}$.
/// When $A$ is non-singular and the Petrov–Galerkin conditions hold without
/// breakdown ($\rho_k \neq 0$, $\omega_k \neq 0$), the algorithm is
/// mathematically equivalent to coupled Lanczos with stabilisation, and
/// the residual norm converges at least as fast as Bi-CG.
///
/// **Proof sketch**: The biorthogonality relation
/// $\langle \tilde{r}_0, r_k \rangle = \rho_k$ and the stabilisation step
/// $r_k \leftarrow s - \omega_k A s$ with $\omega_k = (As)^T s / \|As\|^2$
/// minimise $\|r_k\|_2$ over the one-dimensional affine subspace at each
/// half-step, preventing the erratic convergence of plain Bi-CG.
///
/// **Reference**: Van der Vorst (1992), Theorem 3.1.
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

fn is_finite_scalar<T: RealField + Copy + NumericElement>(x: T) -> bool {
    NumericElement::to_f64(x).is_finite()
}

impl<T: RealField + Copy + NumericElement> BiCGSTAB<T> {
    /// Create new BiCGSTAB solver
    pub const fn new(config: IterativeSolverConfig<T>) -> Self {
        Self { config }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self
    where
        T: FloatElement,
    {
        Self::new(IterativeSolverConfig::default())
    }

    /// Solve with left preconditioning
    pub fn solve_preconditioned<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &Array1<T>,
        preconditioner: &P,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        let n = vector_len(b);
        validate_vector_len("BiCGSTAB solution", x, n)?;
        let a_size = a.size();
        if a_size != 0 && a_size != n {
            return Err(Error::InvalidConfiguration(format!(
                "Operator size ({a_size}) doesn't match RHS vector ({n})"
            )));
        }

        let mut r = Array1::zeros([n]);
        let mut r0_hat = Array1::zeros([n]);
        let mut p = Array1::zeros([n]);
        let mut v = Array1::zeros([n]);
        let mut s = Array1::zeros([n]);
        let mut t = Array1::zeros([n]);
        let mut z = Array1::zeros([n]);
        let mut z2 = Array1::zeros([n]);
        let mut ax = Array1::zeros([n]);

        a.apply(x, &mut ax)?;
        assign_residual(&mut r, b, &ax);

        let initial_residual_norm = norm(&r);
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        if self.is_converged(initial_residual_norm) {
            tracing::debug!("BiCGSTAB converged at initial guess");
            return Ok(monitor);
        }

        let epsilon = <T as RealField>::EPSILON;
        let breakdown_tolerance = epsilon * epsilon;

        copy_array(&r, &mut r0_hat);

        let mut rho = <T as NumericElement>::ONE;
        let mut alpha = <T as NumericElement>::ONE;
        let mut omega = <T as NumericElement>::ONE;

        for _iter in 0..self.config.max_iterations {
            let rho_new = dot(&r0_hat, &r);
            let rho_scale = norm(&r0_hat) * norm(&r);

            if NumericElement::abs(rho_new)
                < breakdown_tolerance * (<T as NumericElement>::ONE + rho_scale)
            {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }

            if NumericElement::abs(omega)
                < breakdown_tolerance * (<T as NumericElement>::ONE + NumericElement::abs(alpha))
            {
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: NumericElement::to_f64(norm(&r)),
                    },
                ));
            }

            let beta = (rho_new / rho) * (alpha / omega);

            axpy(&mut p, -omega, &v);
            scale_add(&mut p, beta, &r);

            preconditioner.apply_to(&p, &mut z)?;
            a.apply(&z, &mut v)?;

            let r0_hat_dot_v = dot(&r0_hat, &v);
            let r0_hat_dot_v_scale = norm(&r0_hat) * norm(&v);
            if NumericElement::abs(r0_hat_dot_v)
                < breakdown_tolerance * (<T as NumericElement>::ONE + r0_hat_dot_v_scale)
            {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            alpha = rho_new / r0_hat_dot_v;
            if !is_finite_scalar(alpha) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }

            copy_array(&r, &mut s);
            axpy(&mut s, -alpha, &v);

            let s_norm = norm(&s);
            if self.is_converged(s_norm) {
                axpy(x, alpha, &z);
                monitor.record_residual(s_norm);
                return Ok(monitor);
            }

            preconditioner.apply_to(&s, &mut z2)?;
            a.apply(&z2, &mut t)?;

            let t_dot_t = dot(&t, &t);
            let t_dot_t_scale = norm(&t) * norm(&t);
            if NumericElement::abs(t_dot_t)
                < breakdown_tolerance * (<T as NumericElement>::ONE + t_dot_t_scale)
            {
                axpy(x, alpha, &z);
                a.apply(x, &mut ax)?;
                assign_residual(&mut r, b, &ax);
                let final_norm = norm(&r);
                monitor.record_residual(final_norm);
                if self.is_converged(final_norm) {
                    return Ok(monitor);
                }
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: NumericElement::to_f64(final_norm),
                    },
                ));
            }

            omega = dot(&t, &s) / t_dot_t;
            if !is_finite_scalar(omega) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }
            axpy(x, alpha, &z);
            axpy(x, omega, &z2);

            copy_array(&s, &mut r);
            axpy(&mut r, -omega, &t);

            let residual_norm = norm(&r);
            if !is_finite_scalar(residual_norm) {
                return Err(Error::Convergence(ConvergenceErrorKind::InvalidValue));
            }
            monitor.record_residual(residual_norm);

            if self.is_converged(residual_norm) {
                return Ok(monitor);
            }

            if NumericElement::abs(omega)
                < breakdown_tolerance * (<T as NumericElement>::ONE + residual_norm)
            {
                return Err(Error::Convergence(
                    ConvergenceErrorKind::StagnatedResidual {
                        residual: NumericElement::to_f64(residual_norm),
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
        b: &Array1<T>,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use super::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }

    fn is_converged(&self, residual_norm: T) -> bool {
        residual_norm < self.config.tolerance
    }
}

impl<T: RealField + Debug + Copy + NumericElement> Configurable<T> for BiCGSTAB<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + NumericElement> IterativeLinearSolver<T> for BiCGSTAB<T> {
    fn solve<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &Array1<T>,
        x: &mut Array1<T>,
        preconditioner: Option<&P>,
    ) -> Result<ConvergenceMonitor<T>> {
        if let Some(p) = preconditioner {
            self.solve_preconditioned(a, b, p, x)
        } else {
            self.solve_unpreconditioned(a, b, x)
        }
    }
}

impl<T: RealField + Copy + FloatElement + NumericElement> super::traits::LinearSolver<T>
    for BiCGSTAB<T>
{
    fn solve_system(
        &self,
        a: &dyn LinearOperator<T>,
        b: &Array1<T>,
        x0: Option<&Array1<T>>,
    ) -> Result<Array1<T>> {
        let mut x = if let Some(initial) = x0 {
            initial.clone()
        } else {
            Array1::zeros(b.shape())
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
    use leto_ops::CsrMatrix;

    fn array(values: Vec<f64>) -> Array1<f64> {
        Array1::from_shape_vec([values.len()], values).expect("valid Leto vector shape")
    }

    fn assert_solves(a: &CsrMatrix<f64>, x: &Array1<f64>, b: &Array1<f64>, epsilon: f64) {
        let mut ax = Array1::zeros([vector_len(b)]);
        a.apply(x, &mut ax).expect("operator application");
        for idx in 0..vector_len(b) {
            assert_relative_eq!(ax[idx], b[idx], epsilon = epsilon);
        }
    }

    fn assert_invalid_configuration<T>(result: Result<T>, expected_message: &str) {
        match result {
            Err(Error::InvalidConfiguration(message)) => assert_eq!(message, expected_message),
            Err(error) => panic!("expected invalid configuration, got {error:?}"),
            Ok(_) => panic!("expected invalid configuration error"),
        }
    }

    fn assert_max_iterations<T>(result: Result<T>, expected_max: usize) {
        match result {
            Err(Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { max })) => {
                assert_eq!(max, expected_max);
            }
            Err(error) => panic!("expected max-iteration convergence error, got {error:?}"),
            Ok(_) => panic!("expected max-iteration convergence error"),
        }
    }

    fn create_nonsymmetric_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![5.0, 1.0, 2.0, 4.0, 1.0, 2.0, 3.0];

        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
    }

    fn create_diagonal_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 1, 2, 3];
        let col_indices = vec![0, 1, 2];
        let values = vec![2.0, 3.0, 4.0];

        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
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
        let b = array(vec![6.0, 11.0, 8.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(
            result.is_ok(),
            "BiCGSTAB should solve nonsymmetric CSR system"
        );
        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_nonsymmetric_matrix();
        let b = array(vec![6.0, 11.0, 8.0]);
        let mut x = array(vec![0.1, 0.2, 0.3]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(
            result.is_ok(),
            "BiCGSTAB should converge from nonzero initial guess"
        );
        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let a = create_diagonal_matrix();
        let b = array(vec![2.0, 3.0, 4.0]);
        let mut x = Array1::zeros([3]);
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
        let b = array(vec![1.0, 2.0]);
        let mut x = Array1::zeros([2]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert_invalid_configuration(result, "Operator size (3) doesn't match RHS vector (2)");
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_nonsymmetric_matrix();
        let b = array(vec![6.0, 11.0, 8.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(
            result.is_ok(),
            "BiCGSTAB should converge with tight tolerance"
        );
        assert_solves(&a, &x, &b, 1e-12);
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_nonsymmetric_matrix();
        let b = array(vec![6.0, 11.0, 8.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::<f64> {
            max_iterations: 1,
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert_max_iterations(result, 1);
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
        let b = Array1::from_shape_vec([3], vec![6.0, 11.0, 8.0]).unwrap();
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);

        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok(), "linear-solver facade should converge");

        let x = result.unwrap();
        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_already_converged() {
        let a = create_diagonal_matrix();
        let b = array(vec![2.0, 3.0, 4.0]);
        let mut x = array(vec![1.0, 1.0, 1.0]);
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
