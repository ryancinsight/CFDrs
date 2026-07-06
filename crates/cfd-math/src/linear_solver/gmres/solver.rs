//! GMRES solver implementation with Arnoldi iteration and Givens rotations

use super::super::config::IterativeSolverConfig;
use super::super::traits::{
    Configurable, ConvergenceMonitor, IterativeLinearSolver, LinearOperator, Preconditioner,
};
use super::{arnoldi, givens};
use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};
use std::fmt::Debug;
use std::sync::Mutex;

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
/// # Theorem — GMRES Krylov Optimality (Saad & Schultz 1986)
///
/// GMRES finds the iterate $x_m \in x_0 + K_m(A, r_0)$ that minimises the
/// 2-norm of the residual over the $m$-th Krylov subspace:
///
/// ```text
/// x_m = x_0 + argmin_{y ∈ K_m(A, r_0)} ‖r_0 − A y‖_2
/// ```
///
/// The residual satisfies the polynomial bound:
///
/// ```text
/// ‖r_m‖ / ‖r_0‖ ≤ inf_{p ∈ Π_m, p(0)=1} max_{λ ∈ σ(A)} |p(λ)|
/// ```
///
/// where $\sigma(A)$ is the spectrum of $A$. For normal matrices the field of
/// values $W(A) = \text{conv}(\sigma(A))$, tightening the bound. For SPD matrices
/// convergence is guaranteed in at most $n$ steps.
///
/// **Proof sketch.** The Arnoldi process builds an orthonormal basis
/// $V_m$ of $K_m(A, r_0)$, satisfying $A V_m = V_{m+1} \bar{H}_m$ where
/// $\bar{H}_m$ is an $(m+1) \times m$ upper Hessenberg matrix. The least-squares
/// minimisation reduces to $\min_y \|\beta e_1 - \bar{H}_m y\|_2$, solved via
/// Givens rotations in $O(m^2)$ work. The polynomial bound follows because
/// the residual polynomial $p_m(A) r_0$ with $p_m(0) = 1$ is the optimal
/// approximation over $\sigma(A)$ in Chebyshev sense (Saad 2003, Thm 6.10).
///
/// ## References
///
/// - Saad, Y. & Schultz, M. H. (1986). "GMRES: A generalized minimal residual
///   algorithm for solving nonsymmetric linear systems." *SIAM J. Sci. Stat.
///   Comput.* 7(3):856–869.
/// - Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.).
///   SIAM. Chapter 6.
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
    /// Cached workspace to prevent reallocation
    workspace: Mutex<Option<GMRESWorkspace<T>>>,
}

#[derive(Clone)]
struct GMRESWorkspace<T: RealField + Copy> {
    v: Array2<T>,
    h: Array2<T>,
    g: Array1<T>,
    c: Array1<T>,
    s: Array1<T>,
    basis_work: Array1<T>,
    work: Array1<T>,
    precond_work: Array1<T>,
    ax: Array1<T>,
}

#[inline]
fn array_vector_fill<T: Copy>(vector: &mut Array1<T>, value: T) {
    for row in 0..vector.shape()[0] {
        vector[row] = value;
    }
}

#[inline]
fn array_matrix_fill<T: Copy>(matrix: &mut Array2<T>, value: T) {
    let [rows, cols] = matrix.shape();
    for row in 0..rows {
        for col in 0..cols {
            matrix[[row, col]] = value;
        }
    }
}

#[inline]
fn array_vector_norm<T: NumericElement>(vector: &Array1<T>) -> T {
    let mut sum = T::ZERO;
    for row in 0..vector.shape()[0] {
        sum += vector[row] * vector[row];
    }
    sum.sqrt()
}

#[inline]
fn array_vector_sub_assign<T: NumericElement>(target: &mut Array1<T>, rhs: &Array1<T>) {
    for row in 0..target.shape()[0] {
        target[row] -= rhs[row];
    }
}

impl<T: RealField + Copy + FloatElement + Debug> GMRES<T> {
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
        b: &Array1<T>,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        use crate::linear_solver::preconditioners::IdentityPreconditioner;
        let preconditioner = IdentityPreconditioner;
        self.solve_preconditioned(a, b, &preconditioner, x)
    }

    /// Solve with left preconditioning using GMRES(m) algorithm
    pub fn solve_preconditioned<Op: LinearOperator<T> + ?Sized, P: Preconditioner<T>>(
        &self,
        a: &Op,
        b: &Array1<T>,
        preconditioner: &P,
        x: &mut Array1<T>,
    ) -> Result<ConvergenceMonitor<T>> {
        let n = b.shape()[0];
        let a_size = a.size();
        if a_size != 0 && a_size != n {
            return Err(Error::InvalidConfiguration(format!(
                "Operator size ({a_size}) doesn't match RHS vector ({n})"
            )));
        }

        let m = self.restart_dim;

        let mut cache_ref = self.workspace.lock().unwrap();
        if cache_ref
            .as_ref()
            .is_none_or(|cache| cache.v.shape() != [n, m + 1])
        {
            *cache_ref = Some(GMRESWorkspace {
                v: Array2::zeros([n, m + 1]),
                h: Array2::zeros([m + 1, m]),
                g: Array1::zeros([m + 1]),
                c: Array1::zeros([m]),
                s: Array1::zeros([m]),
                basis_work: Array1::zeros([n]),
                work: Array1::zeros([n]),
                precond_work: Array1::zeros([n]),
                ax: Array1::zeros([n]),
            });
        }
        let ws = cache_ref.as_mut().unwrap();

        // 1. Initial residual: r0 = b - A*x
        a.apply(x, &mut ws.ax)?;
        let mut r0 = b.clone();
        array_vector_sub_assign(&mut r0, &ws.ax);

        // Apply preconditioning to initial residual if needed (Left Preconditioning)
        preconditioner.apply_to(&r0, &mut ws.work)?;
        let beta = array_vector_norm(&ws.work);

        let r0_norm = array_vector_norm(&r0);
        if r0_norm < self.config.tolerance {
            return Ok(ConvergenceMonitor::new(r0_norm));
        }

        if beta <= <T as RealField>::EPSILON {
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
                let mut r_restart = b.clone();
                array_vector_sub_assign(&mut r_restart, &ws.ax);
                preconditioner.apply_to(&r_restart, &mut ws.work)?;
                array_vector_norm(&ws.work)
            };

            if beta_restart <= <T as RealField>::EPSILON {
                return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
            }

            let inv_beta = <T as NumericElement>::ONE / beta_restart;
            for row in 0..n {
                ws.v[[row, 0]] = ws.work[row] * inv_beta;
            }

            array_matrix_fill(&mut ws.h, <T as NumericElement>::ZERO);
            array_vector_fill(&mut ws.g, <T as NumericElement>::ZERO);
            array_vector_fill(&mut ws.c, <T as NumericElement>::ZERO);
            array_vector_fill(&mut ws.s, <T as NumericElement>::ZERO);
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
                    &mut ws.basis_work,
                    &mut ws.work,
                    Some(preconditioner),
                    Some(&mut ws.precond_work),
                )?;
                iterations_used += 1;

                // Apply previous Givens rotations to new column of H
                givens::apply_previous_rotations(&mut ws.h, &ws.c, &ws.s, k);

                // Compute new Givens rotation to zero out H(k+1, k)
                let (ck, sk) = givens::compute_rotation(ws.h[[k, k]], ws.h[[k + 1, k]]);
                ws.c[k] = ck;
                ws.s[k] = sk;

                // Apply new Givens rotation to H and g
                givens::apply_new_rotation(&mut ws.h, &mut ws.g, ck, sk, k);

                // Check convergence using residual norm estimate
                let residual_estimate = NumericElement::abs(ws.g[k + 1]);
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
                for row in 0..n {
                    x[row] += y[i] * ws.v[[row, i]];
                }
            }

            a.apply(x, &mut ws.ax)?;
            let mut r_check = b.clone();
            array_vector_sub_assign(&mut r_check, &ws.ax);
            if array_vector_norm(&r_check) < self.config.tolerance {
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

impl<T: RealField + Copy + FloatElement + Debug> Configurable<T> for GMRES<T> {
    type Config = IterativeSolverConfig<T>;

    fn config(&self) -> &Self::Config {
        &self.config
    }
}

impl<T: RealField + Debug + Copy + FloatElement> IterativeLinearSolver<T> for GMRES<T> {
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

impl<T: RealField + Copy + FloatElement + Debug> super::super::traits::LinearSolver<T>
    for GMRES<T>
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
            None::<&crate::linear_solver::preconditioners::IdentityPreconditioner>,
        )?;
        Ok(x)
    }
}
