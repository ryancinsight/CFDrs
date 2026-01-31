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
use std::sync::Mutex;

/// Workspace for GMRES solver to avoid repeated allocations
struct GmresWorkspace<T: RealField + Copy> {
    v: DMatrix<T>,
    h: DMatrix<T>,
    g: DVector<T>,
    c: DVector<T>,
    s: DVector<T>,
    work: DVector<T>,
    precond_work: DVector<T>,
    ax: DVector<T>,
    v_col: DVector<T>,
    r: DVector<T>,
    n: usize,
    m: usize,
}

impl<T: RealField + Copy> GmresWorkspace<T> {
    fn new(n: usize, m: usize) -> Self {
        Self {
            v: DMatrix::zeros(n, m + 1),
            h: DMatrix::zeros(m + 1, m),
            g: DVector::zeros(m + 1),
            c: DVector::zeros(m),
            s: DVector::zeros(m),
            work: DVector::zeros(n),
            precond_work: DVector::zeros(n),
            ax: DVector::zeros(n),
            v_col: DVector::zeros(n),
            r: DVector::zeros(n),
            n,
            m,
        }
    }
}

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
    /// Reusable workspace to avoid allocations
    workspace: Mutex<Option<GmresWorkspace<T>>>,
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
            workspace: Mutex::new(None),
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
            return Err(Error::InvalidConfiguration(format!(
                "Operator size ({a_size}) doesn't match RHS vector ({n})"
            )));
        }

        let m = self.restart_dim;

        // Workspace vectors and matrices
        let mut workspace_lock = self.workspace.lock().unwrap();
        if workspace_lock.as_ref().map_or(true, |ws| ws.n != n || ws.m != m) {
            *workspace_lock = Some(GmresWorkspace::new(n, m));
        }
        let ws = workspace_lock.as_mut().unwrap();

        // 1. Initial residual: r0 = b - A*x
        a.apply(x, &mut ws.ax)?;
        ws.r.copy_from(b);
        ws.r -= &ws.ax;

        // Apply preconditioning to initial residual if needed (Left Preconditioning)
        preconditioner.apply_to(&ws.r, &mut ws.work)?;
        let beta = ws.work.norm();

        let r0_norm = ws.r.norm();
        if r0_norm < self.config.tolerance {
            return Ok(ConvergenceMonitor::new(r0_norm));
        }

        if beta <= T::default_epsilon() {
            return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
        }

        let mut monitor = ConvergenceMonitor::new(beta);

        let mut iterations_used = 0usize;
        let mut is_first_restart = true;
        while iterations_used < self.config.max_iterations {
            let beta_restart = if is_first_restart {
                beta
            } else {
                a.apply(x, &mut ws.ax)?;
                ws.r.copy_from(b);
                ws.r -= &ws.ax;
                preconditioner.apply_to(&ws.r, &mut ws.work)?;
                ws.work.norm()
            };

            if beta_restart <= T::default_epsilon() {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }

            {
                let mut v0 = ws.v.column_mut(0);
                v0.copy_from(&ws.work);
                v0 *= T::one() / beta_restart;
            }

            ws.h.fill(T::zero());
            ws.g.fill(T::zero());
            ws.c.fill(T::zero());
            ws.s.fill(T::zero());
            ws.g[0] = beta_restart;

            // 3. Arnoldi iterations
            let mut converged_iter = None;
            let remaining = self.config.max_iterations - iterations_used;
            let inner_iters = m.min(remaining);
            for k in 0..inner_iters {
                // Arnoldi step: build orthonormal basis
                arnoldi::arnoldi_iteration(
                    a,
                    &mut ws.v,
                    &mut ws.h,
                    k,
                    &mut ws.work,
                    Some(preconditioner),
                    Some(&mut ws.precond_work),
                    &mut ws.v_col,
                )?;
                iterations_used += 1;

                // Apply previous Givens rotations to new column of H
                givens::apply_previous_rotations(&mut ws.h, &ws.c, &ws.s, k);

                // Compute new Givens rotation to zero out H(k+1, k)
                let (ck, sk) = givens::compute_rotation(ws.h[(k, k)], ws.h[(k + 1, k)]);
                ws.c[k] = ck;
                ws.s[k] = sk;

                // Apply new Givens rotation to H and g
                givens::apply_new_rotation(&mut ws.h, &mut ws.g, ck, sk, k);

                // Check convergence using residual norm estimate
                let residual_estimate = ws.g[k + 1].abs();
                monitor.record_residual(residual_estimate);

                if residual_estimate < self.config.tolerance {
                    converged_iter = Some(k + 1);
                    break;
                }
            }

            // 4. Update solution: x = x + V_k * y_k
            let k_final = converged_iter.unwrap_or(inner_iters);
            if k_final == 0 {
                break;
            }
            let y = givens::solve_upper_triangular(&ws.h, &ws.g, k_final)?;

            for i in 0..k_final {
                x.axpy(y[i], &ws.v.column(i), T::one());
            }

            a.apply(x, &mut ws.ax)?;
            ws.r.copy_from(b);
            ws.r -= &ws.ax;
            if ws.r.norm() < self.config.tolerance {
                return Ok(monitor);
            }

            is_first_restart = false;
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.max_iterations,
            },
        ))
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
