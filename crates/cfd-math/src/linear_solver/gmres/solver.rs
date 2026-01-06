//! GMRES solver implementation with Arnoldi iteration and Givens rotations

use super::super::config::IterativeSolverConfig;
use super::super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use super::{arnoldi, givens};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DMatrix, DVector, RealField};
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
/// # Convergence Theory
///
/// ## Field of Values Convergence Bound (Saad & Schultz, 1986)
///
/// GMRES convergence is governed by the field of values (numerical range) of the matrix A.
/// The residual norm satisfies:
/// ||r_k|| ≤ κ(V_k) * inf_{p∈Π_k} max_{z∈W(A)} |p(z)| / min_{z∈W(A)} |p(z)|
///
/// where:
/// - Π_k is the set of polynomials of degree ≤ k
/// - W(A) is the field of values of A
/// - κ(V_k) is the condition number of the Vandermonde matrix
///
/// **Theorem (Saad & Schultz, 1986)**: For any matrix A, the GMRES residual satisfies:
/// ||r_m|| / ||r0|| ≤ inf_{p∈Π_m, p(0)=1} max_{λ∈σ(A)} |p(λ)|
///
/// where σ(A) is the spectrum of A. For normal matrices, W(A) = σ(A).
///
/// For symmetric positive definite matrices, convergence is guaranteed in at most n steps.
///
/// ## Optimal Polynomial Approximation
///
/// GMRES finds the vector in K_m(A,r0) that minimizes the residual norm:
/// x_m = x0 + argmin_{y∈K_m} ||r0 - A*y||
///
/// This corresponds to finding the polynomial p_m(z) = 1 - z * q_{m-1}(z) where q_{m-1}
/// minimizes the maximum of |p_m(z)| over the field of values W(A).
///
/// The convergence factor satisfies:
/// ||r_m|| / ||r0|| ≤ inf_{p∈Π_m, p(0)=1} max_{z∈W(A)} |p(z)|
///
/// ## Restart Parameter Justification
///
/// The restart dimension m=30 is chosen as a practical compromise:
/// - Memory usage: O(n*m) for basis vectors, O(m²) for Hessenberg matrix
/// - For CFD applications, m=20-50 typically provides good convergence
/// - Larger m reduces restart overhead but increases memory usage
/// - m=30 balances computational efficiency with convergence speed
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
        Self {
            config,
            restart_dim,
        }
    }

    /// Create with default configuration
    #[must_use]
    pub fn default() -> Self {
        Self::new(IterativeSolverConfig::default(), 30)
    }

    /// Solve without preconditioning using GMRES(m) algorithm
    pub fn solve_unpreconditioned<Op: LinearOperator<T> + ?Sized>(
        &self,
        a: &Op,
        b: &DVector<T>,
        x: &mut DVector<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use crate::linear_solver::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }

    /// Solve with left preconditioning using GMRES(m) algorithm
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

        let m = self.restart_dim;

        // Workspace vectors and matrices
        let mut v = DMatrix::zeros(n, m + 1);
        let mut h = DMatrix::zeros(m + 1, m);
        let mut g = DVector::zeros(m + 1);
        let mut c = DVector::zeros(m);
        let mut s = DVector::zeros(m);
        let mut work = DVector::zeros(n);
        let mut precond_work = DVector::zeros(n);
        let mut ax = DVector::zeros(n);

        // 1. Initial residual: r0 = b - A*x
        a.apply(x, &mut ax)?;
        let mut r0 = b.clone();
        r0 -= &ax;

        // Apply preconditioning to initial residual if needed (Left Preconditioning)
        preconditioner.apply_to(&r0, &mut work)?;
        let beta = work.norm();

        let mut monitor = ConvergenceMonitor::new(beta);

        if beta < self.config.tolerance {
            return Ok(monitor);
        }

        for _outer_iter in 0..=(self.config.max_iterations / m) {
            if _outer_iter > 0 {
                // Compute residual for restart
                a.apply(x, &mut ax)?;
                let mut r_restart = b.clone();
                r_restart -= &ax;
                preconditioner.apply_to(&r_restart, &mut work)?;
                let beta_restart = work.norm();
                monitor.record_residual(beta_restart);

                if beta_restart < self.config.tolerance {
                    return Ok(monitor);
                }
                
                // Initialize first basis vector for restart
                v.column_mut(0).copy_from(&(work.clone() / beta_restart));
                g.fill(T::zero());
                g[0] = beta_restart;
            } else {
                // 2. Initialize first basis vector: v1 = r0 / beta
                v.column_mut(0).copy_from(&(work.clone() / beta));
                g.fill(T::zero());
                g[0] = beta;
            }

            // 3. Arnoldi iterations
            let mut converged_iter = None;
            for k in 0..m {
                // Arnoldi step: build orthonormal basis
                arnoldi::arnoldi_iteration(
                    a,
                    &mut v,
                    &mut h,
                    k,
                    &mut work,
                    Some(preconditioner),
                    Some(&mut precond_work),
                )?;

                // Apply previous Givens rotations to new column of H
                givens::apply_previous_rotations(&mut h, &c, &s, k);

                // Compute new Givens rotation to zero out H(k+1, k)
                let (ck, sk) = givens::compute_rotation(h[(k, k)], h[(k + 1, k)]);
                c[k] = ck;
                s[k] = sk;

                // Apply new Givens rotation to H and g
                givens::apply_new_rotation(&mut h, &mut g, ck, sk, k);

                // Check convergence using residual norm estimate
                let residual_estimate = g[k + 1].abs();
                monitor.record_residual(residual_estimate);

                if residual_estimate < self.config.tolerance {
                    converged_iter = Some(k + 1);
                    break;
                }
            }

            // 4. Update solution: x = x + V_k * y_k
            let k_final = converged_iter.unwrap_or(m);
            let y = givens::solve_upper_triangular(&h, &g, k_final)?;

            for i in 0..k_final {
                x.axpy(y[i], &v.column(i), T::one());
            }

            if converged_iter.is_some() {
                return Ok(monitor);
            }
        }

        Err(Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded {
            max: self.config.max_iterations,
        }))
    }
}

impl<T: RealField + Copy + FromPrimitive + Debug> Configurable<T> for GMRES<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + FromPrimitive> IterativeLinearSolver<T> for GMRES<T> {
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

impl<T: RealField + Copy + FromPrimitive + Debug> super::super::traits::LinearSolver<T>
    for GMRES<T>
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
            None::<&crate::linear_solver::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}

#[cfg(test)]
mod tests {
    use super::super::super::preconditioners::IdentityPreconditioner;
    use super::*;
    use nalgebra_sparse::{CooMatrix, CsrMatrix};

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
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        // Expected solution: x = [1, 2, 3, 4]
        let expected = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let error = (&x - &expected).norm();
        assert!(error < 1e-9, "Solution error: {error}");
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
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        // Verify Ax = b
        let ax = &a * &x;
        let residual = (&ax - &b).norm();
        assert!(residual < 1e-9, "Residual: {residual}");
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
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        assert!(x.norm() < 1e-14, "Solution should be zero for zero RHS");
    }

    #[test]
    fn test_gmres_restart_dimension() {
        let n = 5;
        let mut coo = CooMatrix::new(n, n);
        // Tridiagonal matrix
        for i in 0..n {
            coo.push(i, i, 4.0);
            if i > 0 {
                coo.push(i, i - 1, 1.0);
            }
            if i < n - 1 {
                coo.push(i, i + 1, 1.0);
            }
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 3); // Small restart dimension
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        // Verify solution
        let ax = &a * &x;
        let residual = (&ax - &b).norm();
        assert!(residual < 1e-9, "Residual: {residual}");
    }

    #[test]
    fn test_gmres_with_initial_guess() {
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
        let mut x = DVector::from_vec(vec![0.5, 1.5, 2.5, 3.5]); // Good initial guess

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        // Expected solution: x = [1, 2, 3, 4]
        let expected = DVector::from_vec(vec![1.0, 2.0, 3.0, 4.0]);
        let error = (&x - &expected).norm();
        assert!(error < 1e-9, "Solution error: {error}");
    }

    #[test]
    fn test_gmres_larger_nonsymmetric() {
        // 5x5 nonsymmetric matrix
        let n = 5;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, 5.0);
            if i > 0 {
                coo.push(i, i - 1, 2.0);
            }
            if i < n - 1 {
                coo.push(i, i + 1, 1.0);
            }
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![6.0, 11.0, 11.0, 11.0, 8.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner;

        solver
            .solve_preconditioned(&a, &b, &precond, &mut x)
            .unwrap();

        // Verify solution
        let ax = &a * &x;
        let residual = (&ax - &b).norm();
        assert!(residual < 1e-9, "Residual: {residual}");
    }

    #[test]
    fn test_gmres_convergence_tight_tolerance() {
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-14).with_max_iterations(100); // Very tight tolerance
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_ok());
    }

    #[test]
    fn test_gmres_max_iterations_exceeded() {
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 4.0, 9.0, 16.0]);
        let mut x = DVector::zeros(n);

        let config = IterativeSolverConfig::new(1e-14).with_max_iterations(1); // Too few iterations
        let solver = GMRES::new(config, 1);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err()); // Should fail to converge
    }

    #[test]
    fn test_gmres_dimension_mismatch() {
        let n = 4;
        let mut coo = CooMatrix::new(n, n);
        for i in 0..n {
            coo.push(i, i, (i + 1) as f64);
        }
        let a = CsrMatrix::from(&coo);

        let b = DVector::from_vec(vec![1.0, 4.0]); // Wrong size!
        let mut x = DVector::zeros(2);

        let config = IterativeSolverConfig::new(1e-10).with_max_iterations(100);
        let solver = GMRES::new(config, 10);
        let precond = IdentityPreconditioner;

        let result = solver.solve_preconditioned(&a, &b, &precond, &mut x);
        assert!(result.is_err());
    }

    #[test]
    fn test_gmres_configurable_trait() {
        let config = IterativeSolverConfig::new(1e-8).with_max_iterations(200);
        let solver = GMRES::<f64>::new(config, 20);

        let retrieved = solver.config();
        assert!((retrieved.tolerance - 1e-8).abs() < 1e-10);
        assert_eq!(retrieved.max_iterations, 200);
    }

    #[test]
    #[should_panic(expected = "GMRES restart dimension must be positive")]
    fn test_gmres_zero_restart_dim_panics() {
        let config = IterativeSolverConfig::default();
        let _solver = GMRES::<f64>::new(config, 0); // Should panic
    }
}
