//! Preconditioned Conjugate Gradient (PCG) solver implementation.
//!
//! ## Literature References
//!
//! - Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients for solving linear equations.
//! - Saad, Y. (2003): *Iterative Methods for Sparse Linear Systems*, SIAM. Section 6.7.
//! - Barrett et al. (1994): *Templates for the Solution of Linear Systems*, SIAM.
//! - Golub, G. H., & Van Loan, C. F. (2013): *Matrix Computations* (4th ed.).
//! - Meurant, G. (1999): *Computer Solution of Large Linear Systems*.

use super::array_ops::{
    assign_residual, axpy, copy_array, dot, norm, scale_add, validate_vector_len, vector_len,
};
use super::config::IterativeSolverConfig;
use super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;
use std::fmt::Debug;
use std::sync::Mutex;

/// Workspace for CG solver to avoid repeated allocations
struct Workspace<T: RealField + Copy> {
    r: Array1<T>,
    p: Array1<T>,
    z: Array1<T>,
    ap: Array1<T>,
    ax: Array1<T>,
}

/// Preconditioned Conjugate Gradient solver
///
/// # Theorem (CG Convergence in the A-norm)
///
/// For an SPD matrix $A \in \mathbb{R}^{n \times n}$ with eigenvalues
/// $0 < \lambda_1 \le \cdots \le \lambda_n$, the CG iterates satisfy
///
/// $$\|x_k - x^*\|_A \le 2 \left(\frac{\sqrt{\kappa} - 1}{\sqrt{\kappa} + 1}\right)^k \|x_0 - x^*\|_A$$
///
/// where $\kappa = \lambda_n / \lambda_1$ is the spectral condition number.
///
/// **Proof sketch**: CG minimises $\|x_k - x^*\|_A$ over the Krylov subspace
/// $x_0 + \mathcal{K}_k(A, r_0)$. Among all polynomials $p_k$ of degree $k$
/// with $p_k(0) = 1$, the Chebyshev polynomial $T_k$ on $[\lambda_1, \lambda_n]$
/// yields the tightest uniform bound. Evaluating $\max_{\lambda \in [\lambda_1,
/// \lambda_n]} |T_k(\lambda) / T_k(0)|$ gives the stated contraction factor.
///
/// **Reference**: Hestenes & Stiefel (1952); Saad (2003) §6.7.
pub struct ConjugateGradient<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
    /// Thread-safe workspace for reusing allocations
    workspace: Mutex<Option<Workspace<T>>>,
}

impl<T: RealField + Copy + NumericElement> ConjugateGradient<T> {
    /// Create new CG solver
    pub const fn new(config: IterativeSolverConfig<T>) -> Self {
        Self {
            config,
            workspace: Mutex::new(None),
        }
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
        validate_vector_len("CG solution", x, n)?;
        let a_size = a.size();
        if a_size != 0 && a_size != n {
            return Err(Error::InvalidConfiguration(format!(
                "Operator size ({a_size}) doesn't match RHS vector ({n})"
            )));
        }

        // Acquire workspace
        let mut guard = self.workspace.lock().unwrap();

        // Initialize or resize workspace if needed
        if guard.as_ref().is_none_or(|ws| vector_len(&ws.r) != n) {
            *guard = Some(Workspace {
                r: Array1::zeros([n]),
                p: Array1::zeros([n]),
                z: Array1::zeros([n]),
                ap: Array1::zeros([n]),
                ax: Array1::zeros([n]),
            });
        }

        let Workspace { r, p, z, ap, ax } = guard.as_mut().unwrap();

        a.apply(x, ax)?;
        assign_residual(r, b, ax);

        let initial_residual_norm = norm(r);
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        if initial_residual_norm < self.config.tolerance {
            return Ok(monitor);
        }

        preconditioner.apply_to(r, z)?;
        copy_array(z, p);

        let epsilon = <T as RealField>::EPSILON;
        let breakdown_tolerance = epsilon * epsilon;

        let mut r_dot_z = dot(r, z);
        let r_dot_z_scale = norm(r) * norm(z);
        if NumericElement::abs(r_dot_z)
            < breakdown_tolerance * (<T as NumericElement>::ONE + r_dot_z_scale)
        {
            return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
        }
        if r_dot_z < <T as NumericElement>::ZERO {
            return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
        }

        for _iter in 0..self.config.max_iterations {
            a.apply(p, ap)?;

            let p_dot_ap = dot(p, ap);
            let p_dot_ap_scale = norm(p) * norm(ap);
            if NumericElement::abs(p_dot_ap)
                < breakdown_tolerance * (<T as NumericElement>::ONE + p_dot_ap_scale)
            {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            if p_dot_ap < <T as NumericElement>::ZERO {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            let alpha = r_dot_z / p_dot_ap;

            axpy(x, alpha, p);
            axpy(r, -alpha, ap);

            let residual_norm = norm(r);
            monitor.record_residual(residual_norm);

            if residual_norm < self.config.tolerance {
                return Ok(monitor);
            }

            preconditioner.apply_to(r, z)?;

            let r_dot_z_new = dot(r, z);
            let r_dot_z_new_scale = norm(r) * norm(z);
            if NumericElement::abs(r_dot_z_new)
                < breakdown_tolerance * (<T as NumericElement>::ONE + r_dot_z_new_scale)
            {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }
            if r_dot_z_new < <T as NumericElement>::ZERO {
                return Err(Error::Numerical(NumericalErrorKind::NotPositiveDefinite));
            }

            let beta = r_dot_z_new / r_dot_z;

            scale_add(p, beta, z);

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
        b: &Array1<T>,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use super::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }
}

impl<T: RealField + Debug + Copy + NumericElement> Configurable<T> for ConjugateGradient<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + NumericElement> IterativeLinearSolver<T>
    for ConjugateGradient<T>
{
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
    for ConjugateGradient<T>
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

    fn create_simple_spd_matrix() -> CsrMatrix<f64> {
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0];

        CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3).expect("Valid CSR matrix")
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
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(
            result.is_ok(),
            "CG should converge from nonzero initial guess"
        );
        assert_solves(&a, &x, &b, 1e-6);

        assert_solves(&a, &x, &b, 1e-6);
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = array(vec![0.1, 0.2, 0.3]);
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
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3)
            .expect("Valid CSR matrix");

        let b = array(vec![5.0, 10.0, 15.0]);
        let mut x = Array1::zeros([3]);
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
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, 3, 3)
            .expect("Valid CSR matrix");

        let b = array(vec![6.0, 9.0, 12.0]);
        let mut x = Array1::zeros([3]);
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
        let b = array(vec![1.0, 2.0]);
        let mut x = Array1::zeros([2]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert_invalid_configuration(result, "Operator size (3) doesn't match RHS vector (2)");
    }

    #[test]
    fn test_convergence_with_tight_tolerance() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::<f64> {
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok(), "CG should converge with tight tolerance");
        assert_solves(&a, &x, &b, 1e-12);
    }

    #[test]
    fn test_max_iterations_exceeded() {
        let a = create_simple_spd_matrix();
        let b = array(vec![1.0, 2.0, 3.0]);
        let mut x = Array1::zeros([3]);
        let config = IterativeSolverConfig::<f64> {
            max_iterations: 1,
            tolerance: 1e-12,
            ..Default::default()
        };
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert_max_iterations(result, 1);
    }

    #[test]
    fn test_solve_larger_system() {
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4];
        let values = vec![
            4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0, 1.0, 1.0, 4.0,
        ];
        let a = CsrMatrix::from_parts(values, col_indices, row_offsets, 5, 5)
            .expect("Valid CSR matrix");

        let b = array(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = Array1::zeros([5]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        assert_solves(&a, &x, &b, 1e-6);
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
        let b = array(vec![1.0, 2.0, 3.0]);
        let config = IterativeSolverConfig::default();
        let solver = ConjugateGradient::new(config);

        let result = solver.solve_system(&a, &b, None);
        assert!(result.is_ok());

        let x = result.unwrap();
        assert_solves(&a, &x, &b, 1e-6);
    }
}
