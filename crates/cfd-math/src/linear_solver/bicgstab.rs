//! BiCGSTAB solver implementation

use super::config::IterativeSolverConfig;
use super::traits::{Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner};
use cfd_core::error::{ConvergenceErrorKind, Error, NumericalErrorKind, Result};
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
///   Bi-CG methods. *Numerical Algorithms*, 7(1), 75-109.
///   See Section 2.1 for BiCGSTAB algorithm derivation.
/// - Saad, Y. (2003). *Iterative methods for sparse linear systems* (2nd ed.).
///   SIAM. Section 7.3: BiCG and Variants, including BiCGSTAB.
/// - Sonneveld, P. (1989). CGS, a fast Lanczos-type solver for nonsymmetric linear systems.
///   *SIAM Journal on Scientific and Statistical Computing*, 10(1), 36-52.
///   Original CGS method, precursor to BiCGSTAB.
/// - Fletcher, R. (1976). Conjugate gradient methods for indefinite systems.
///   In *Numerical Analysis* (pp. 73-89). Springer.
///   BiCG foundation and convergence theory.
pub struct BiCGSTAB<T: RealField + Copy> {
    config: IterativeSolverConfig<T>,
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

    /// Solve with preconditioning using efficient memory management
    ///
    /// # Arguments
    /// * `a` - The coefficient operator
    /// * `b` - The right-hand side vector
    /// * `preconditioner` - The preconditioner to use
    /// * `x` - On entry: initial guess; On exit: solution vector
    ///
    /// # Returns
    /// * `Ok(ConvergenceMonitor)` if converged successfully
    /// * `Err(...)` if failed to converge or numerical breakdown occurred
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
            return Err(Error::InvalidConfiguration(
                format!("Operator size ({a_size}) doesn't match RHS vector ({n})"),
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
        let mut v = DVector::zeros(n);
        let mut s = DVector::zeros(n);
        let mut t = DVector::zeros(n);
        let mut z = DVector::zeros(n);
        let mut z2 = DVector::zeros(n);
        let mut ax = DVector::zeros(n);

        // Compute initial residual: r = b - A*x
        a.apply(x, &mut ax)?;
        r.copy_from(b);
        r -= &ax;

        let initial_residual_norm = r.norm();
        let mut monitor = ConvergenceMonitor::new(initial_residual_norm);

        // Check if already converged
        if self.is_converged(initial_residual_norm) {
            tracing::debug!("BiCGSTAB converged at initial guess");
            return Ok(monitor);
        }

        // Use a robust breakdown tolerance based on machine epsilon
        let epsilon = T::default_epsilon();
        let breakdown_tolerance = epsilon * epsilon; // Simple, robust choice

        r0_hat.copy_from(&r); // Shadow residual

        let mut rho = T::one();
        let mut alpha = T::one();
        let mut omega = T::one();

        for _iter in 0..self.config.max_iterations {
            let rho_new = r0_hat.dot(&r);

            // Primary breakdown condition: rho approaching zero
            if rho_new.abs() < breakdown_tolerance {
                return Err(Error::Numerical(NumericalErrorKind::SingularMatrix));
            }

            let beta = (rho_new / rho) * (alpha / omega);

            p.axpy(-omega, &v, T::one());
            p *= beta;
            p += &r;

            preconditioner.apply_to(&p, &mut z)?;
            a.apply(&z, &mut v)?;

            alpha = rho_new / r0_hat.dot(&v);

            s.copy_from(&r);
            s.axpy(-alpha, &v, T::one());

            preconditioner.apply_to(&s, &mut z2)?;
            a.apply(&z2, &mut t)?;

            let t_dot_t = t.dot(&t);
            if t_dot_t.abs() < breakdown_tolerance {
                x.axpy(alpha, &z, T::one());
                a.apply(x, &mut ax)?;
                r.copy_from(b);
                r -= &ax;
                let final_norm = r.norm();
                monitor.record_residual(final_norm);
                if self.is_converged(final_norm) {
                    return Ok(monitor);
                }
                return Err(Error::Convergence(ConvergenceErrorKind::StagnatedResidual {
                    residual: final_norm.to_subset().unwrap_or(0.0),
                }));
            }

            omega = t.dot(&s) / t_dot_t;
            x.axpy(alpha, &z, T::one());
            x.axpy(omega, &z2, T::one());

            r.copy_from(&s);
            r.axpy(-omega, &t, T::one());

            let residual_norm = r.norm();
            monitor.record_residual(residual_norm);

            if self.is_converged(residual_norm) {
                return Ok(monitor);
            }

            if omega.abs() < breakdown_tolerance {
                return Err(Error::Convergence(ConvergenceErrorKind::StagnatedResidual {
                    residual: residual_norm.to_subset().unwrap_or(0.0),
                }));
            }

            rho = rho_new;
        }

        Err(Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
            max: self.config.max_iterations,
        }))
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
    use crate::sparse::spmv;
    use approx::assert_relative_eq;
    use nalgebra_sparse::CsrMatrix;

    fn create_nonsymmetric_matrix() -> CsrMatrix<f64> {
        // Create a 3x3 nonsymmetric matrix (for testing BiCGSTAB)
        // [5, 1, 0]
        // [2, 4, 1]
        // [0, 2, 3]
        let row_offsets = vec![0, 2, 5, 7];
        let col_indices = vec![0, 1, 0, 1, 2, 1, 2];
        let values = vec![5.0, 1.0, 2.0, 4.0, 1.0, 2.0, 3.0];

        CsrMatrix::try_from_csr_data(3, 3, row_offsets, col_indices, values)
            .expect("Valid CSR matrix")
    }

    fn create_diagonal_matrix() -> CsrMatrix<f64> {
        // Diagonal matrix [2, 0, 0; 0, 3, 0; 0, 0, 4]
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
    fn test_solve_nonsymmetric_system() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        // Verify solution by checking A*x ≈ b
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_diagonal_matrix() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![6.0, 9.0, 12.0]);
        let mut x = DVector::zeros(3);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        // Solution should be [3.0, 3.0, 3.0]
        assert_relative_eq!(x[0], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 3.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 3.0, epsilon = 1e-10);
    }

    #[test]
    fn test_solve_with_initial_guess() {
        let a = create_nonsymmetric_matrix();
        let b = DVector::from_vec(vec![6.0, 11.0, 8.0]);
        let mut x = DVector::from_vec(vec![0.5, 0.5, 0.5]);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        // Verify solution
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_mismatched_dimensions_matrix() {
        let a = create_nonsymmetric_matrix(); // 3x3
        let b = DVector::from_vec(vec![1.0, 2.0]); // Wrong size!
        let mut x = DVector::zeros(2);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_mismatched_dimensions_solution() {
        let a = create_nonsymmetric_matrix(); // 3x3
        let b = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let mut x = DVector::zeros(2); // Wrong size!
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
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-12; // Very tight tolerance
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
        let mut config = IterativeSolverConfig::default();
        config.max_iterations = 1; // Too few iterations
        config.tolerance = 1e-12; // Tight tolerance
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        // Should fail to converge
        assert!(result.is_err());
    }

    #[test]
    fn test_already_converged() {
        let a = create_diagonal_matrix();
        let b = DVector::from_vec(vec![2.0, 3.0, 4.0]);
        let mut x = DVector::from_vec(vec![1.0, 1.0, 1.0]); // Exact solution
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        // Solution should remain unchanged
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[1], 1.0, epsilon = 1e-10);
        assert_relative_eq!(x[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_configurable_trait() {
        let mut config = IterativeSolverConfig::default();
        config.tolerance = 1e-8;
        config.max_iterations = 500;

        let solver = BiCGSTAB::new(config);

        // Test getting config
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
        let mut ax = DVector::zeros(3);
        spmv(&a, &x, &mut ax);
        for i in 0..3 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_solve_larger_nonsymmetric_system() {
        // 5x5 nonsymmetric tridiagonal matrix
        let row_offsets = vec![0, 2, 5, 8, 11, 13];
        let col_indices = vec![
            0, 1, // row 0
            0, 1, 2, // row 1
            1, 2, 3, // row 2
            2, 3, 4, // row 3
            3, 4, // row 4
        ];
        let values = vec![
            5.0, 1.0, // row 0
            2.0, 4.0, 1.0, // row 1
            2.0, 4.0, 1.0, // row 2
            2.0, 4.0, 1.0, // row 3
            2.0, 3.0, // row 4
        ];
        let a = CsrMatrix::try_from_csr_data(5, 5, row_offsets, col_indices, values)
            .expect("Valid CSR matrix");

        let b = DVector::from_vec(vec![6.0, 11.0, 11.0, 11.0, 8.0]);
        let mut x = DVector::zeros(5);
        let config = IterativeSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());

        let mut ax = DVector::zeros(5);
        spmv(&a, &x, &mut ax);
        for i in 0..5 {
            assert_relative_eq!(ax[i], b[i], epsilon = 1e-6);
        }
    }
}
