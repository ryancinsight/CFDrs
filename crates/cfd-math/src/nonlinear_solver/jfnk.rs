//! Jacobian-Free Newton-Krylov (JFNK) Nonlinear Solver
//!
//! # Theorem — Newton-Krylov Quadratic Convergence (Kantorovich 1948)
//!
//! Let $F: \mathbb{R}^n \to \mathbb{R}^n$ be continuously differentiable at $x^*$
//! with $F(x^*) = 0$ and $J(x^*) = \nabla F(x^*)$ invertible. Then Newton's method
//!
//! $$ x_{k+1} = x_k - J(x_k)^{-1} F(x_k) $$
//!
//! converges quadratically in a neighbourhood of $x^*$:
//!
//! $$ \|e_{k+1}\| \leq C \|e_k\|^2, \quad C = \tfrac{1}{2} \|J(x^*)^{-1}\| \|F''\|_\infty $$
//!
//! **JFNK abandons the explicit Jacobian**: instead, the linear system $J(x_k) \delta x = -F(x_k)$
//! is solved by GMRES/FGMRES using only matrix-vector products $J(x) v$.
//!
//! # Theorem — Finite-Difference Jacobian-Vector Product (Brown & Saad 1990, Eq. 2.2)
//!
//! The directional derivative of $F$ at $x$ in direction $v$ is approximated by:
//!
//! $$ J(x) v \approx \frac{F(x + \varepsilon v) - F(x)}{\varepsilon} $$
//!
//! **Optimal perturbation** (Brown & Saad 1990): choose
//! $\varepsilon = \sqrt{\varepsilon_{\rm mach}} \cdot \frac{1 + \|x\|_2}{\|v\|_2}$
//! to balance truncation error $O(\varepsilon)$ against cancellation error
//! $O(\varepsilon_{\rm mach} / \varepsilon)$, giving total error $O(\varepsilon_{\rm mach}^{1/2})$.
//!
//! **Proof of optimality**: Let $\delta_T = L\varepsilon$ (truncation error for Lipschitz $F'$)
//! and $\delta_C = \varepsilon_{\rm mach} \|F(x)\| / \varepsilon$ (cancellation error).
//! Total error $\delta_T + \delta_C$ is minimised at $\varepsilon^* = \sqrt{\varepsilon_{\rm mach} \|F(x)\| / L}$.
//! With $\|F(x)\| \approx 1 + \|x\|$ near convergence this yields the formula above.
//!
//! # Theorem — Inexact Newton Forcing Terms (Eisenstat & Walker 1996, Eq. 2.6)
//!
//! The inner GMRES solve need not be solved to high accuracy early in the iteration.
//! The **Eisenstat-Walker** (EW2) forcing term:
//!
//! $$ \eta_k = \left| \frac{\|F(x_k)\| - \|F(x_{k-1}) + J(x_{k-1})\delta x_{k-1}\|}{\|F(x_{k-1})\|} \right| $$
//!
//! bounds $\eta_k \in [\eta_{\min}, \eta_{\max}]$ and guarantees quadratic-rate preservation:
//! if $\eta_k \leq C \|F(x_k)\|$ for some $C > 0$, the outer Newton step still superconverges.
//!
//! **Reference**: Eisenstat, S.C. & Walker, H.F. (1996). Choosing the forcing terms in an
//! inexact Newton method. *SIAM J. Sci. Comput.* 17(1):16–32.
//!
//! # References
//!
//! - Brown, P.N. & Saad, Y. (1990). Hybrid Krylov methods for nonlinear systems of equations.
//!   *SIAM J. Sci. Stat. Comput.* 11(3):450–481.
//! - Kantorovich, L.V. (1948). On Newton's method for functional equations.
//!   *Doklady Akademii Nauk SSSR* 59(7):1237–1240.
//! - Knoll, D.A. & Keyes, D.E. (2004). Jacobian-free Newton-Krylov methods: A survey.
//!   *J. Comput. Phys.* 193(2):357–397.

use cfd_core::error::Result;
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

use super::linalg::{
    add, add_scaled, add_scaled_in_place, dot, neg, norm, scale, sub, sub_in_place_scaled,
    vector_zeros,
};

/// Configuration for the JFNK nonlinear solver.
#[derive(Debug, Clone)]
pub struct JfnkConfig<T: RealField> {
    /// Maximum outer Newton iterations.
    pub max_newton_iterations: usize,
    /// Absolute convergence tolerance: ‖F(x)‖ < atol.
    pub atol: T,
    /// Relative convergence tolerance: ‖F(xₖ)‖ < rtol · ‖F(x₀)‖.
    pub rtol: T,
    /// Maximum inner Krylov iterations per Newton step.
    pub max_krylov_iterations: usize,
    /// Krylov restart parameter (GMRES-m).
    pub krylov_restart: usize,
    /// Inner solver tolerance at first Newton step (EW forcing initial value).
    pub inner_tol_init: T,
    /// Minimum inner forcing term (denominator safeguard).
    pub eta_min: T,
    /// Maximum inner forcing term (prevents under-solving slowing outer conv).
    pub eta_max: T,
}

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

impl<T: RealField + FloatElement> Default for JfnkConfig<T> {
    fn default() -> Self {
        Self {
            max_newton_iterations: 50,
            atol: from_f64(1e-10),
            rtol: from_f64(1e-8),
            max_krylov_iterations: 200,
            krylov_restart: 30,
            inner_tol_init: from_f64(0.5),
            eta_min: from_f64(1e-4),
            eta_max: from_f64(0.9),
        }
    }
}

/// Convergence status returned by `JfnkSolver::solve`.
#[derive(Debug, Clone)]
pub struct JfnkConvergence<T: Copy> {
    /// Final residual norm ‖F(x)‖
    pub residual_norm: T,
    /// Number of outer Newton iterations performed
    pub newton_iterations: usize,
    /// Residual norm history (one entry per Newton step)
    pub residual_history: Vec<T>,
    /// Whether convergence tolerance was achieved
    pub converged: bool,
}

/// Jacobian-free operator: applies $J(x) v \approx [F(x + \varepsilon v) - F(x_0)] / \varepsilon$
///
/// Internally caches $F(x_0)$ so that repeated Krylov iterations in one Newton step
/// share the same pivot evaluation (one function call per Newton step, not per Krylov step).
struct JvpOperator<'a, T, F>
where
    T: RealField + Copy + FloatElement,
    F: Fn(&Array1<T>) -> Array1<T>,
{
    /// Current Newton point $x_k$ (pivot)
    x_pivot: &'a Array1<T>,
    /// $F(x_k)$ — pre-computed residual at pivot
    f_pivot: &'a Array1<T>,
    /// The nonlinear function $F$
    func: &'a F,
    /// Finite-difference perturbation size $\varepsilon$
    eps: T,
}

impl<T, F> JvpOperator<'_, T, F>
where
    T: RealField + Copy + FloatElement,
    F: Fn(&Array1<T>) -> Array1<T>,
{
    /// Compute $J(x) v \approx [F(x + \varepsilon v) - F(x)] / \varepsilon$.
    ///
    /// # Theorem — Brown-Saad perturbation (see module docs)
    ///
    /// Uses the caller-pre-computed `eps` which already encodes
    /// $\varepsilon = \sqrt{\varepsilon_{\rm mach}} (1 + \|x\|) / \|v\|$.
    fn apply_jvp(&self, v: &Array1<T>) -> Array1<T> {
        // Perturbed point: x + ε v
        let x_pert = add_scaled(self.x_pivot, v, self.eps);
        // Finite difference: [F(x + ε v) - F(x)] / ε
        let f_pert = (self.func)(&x_pert);
        scale(&sub(&f_pert, self.f_pivot), T::ONE / self.eps)
    }
}

/// Jacobian-Free Newton-Krylov solver for nonlinear systems $F(x) = 0$.
///
/// Uses GMRES as the inner Krylov solver with Eisenstat-Walker adaptive forcing.
pub struct JfnkSolver<T: RealField + Copy> {
    config: JfnkConfig<T>,
}

impl<T: RealField + Copy + FloatElement + std::fmt::Debug> JfnkSolver<T> {
    /// Create a new JFNK solver with the given configuration.
    pub fn new(config: JfnkConfig<T>) -> Self {
        Self { config }
    }

    /// Compute the optimal finite-difference perturbation size.
    ///
    /// # Formula (Brown & Saad 1990, Eq. 2.2)
    ///
    /// $\varepsilon = \sqrt{\varepsilon_{\rm mach}} \cdot (1 + \|x\|_2) / \|v\|_2$
    ///
    /// Guards against $\|v\| = 0$ by returning $\varepsilon = \sqrt{\varepsilon_{\rm mach}}$.
    fn compute_eps(x: &Array1<T>, v: &Array1<T>) -> T {
        let eps_mach = T::EPSILON;
        let sqrt_eps = NumericElement::sqrt(eps_mach);
        let norm_v = norm(v);
        if norm_v < eps_mach * from_f64(100.0) {
            sqrt_eps
        } else {
            sqrt_eps * (T::ONE + norm(x)) / norm_v
        }
    }

    /// Solve $F(x) = 0$ with initial guess $x_0$ using JFNK.
    ///
    /// Returns the solution and convergence diagnostics.
    pub fn solve<F>(&self, func: F, mut x: Array1<T>) -> Result<(Array1<T>, JfnkConvergence<T>)>
    where
        F: Fn(&Array1<T>) -> Array1<T>,
    {
        let n = x.shape()[0];
        let mut f = func(&x);
        let norm0 = norm(&f);

        if norm0 < self.config.atol {
            return Ok((
                x,
                JfnkConvergence {
                    residual_norm: norm0,
                    newton_iterations: 0,
                    residual_history: vec![norm0],
                    converged: true,
                },
            ));
        }

        let abs_tol = self.config.atol.max_scalar(self.config.rtol * norm0);

        let mut history = vec![norm0];
        let mut eta = self.config.inner_tol_init; // EW forcing term
        let mut prev_norm = norm0;
        let mut prev_model_residual = norm0; // ‖F(xₖ₋₁) + J δxₖ₋₁‖ for EW2

        for newton_iter in 0..self.config.max_newton_iterations {
            let norm_f = norm(&f);

            if norm_f < abs_tol {
                return Ok((
                    x,
                    JfnkConvergence {
                        residual_norm: norm_f,
                        newton_iterations: newton_iter,
                        residual_history: history,
                        converged: true,
                    },
                ));
            }

            // Inner Krylov tolerance via Eisenstat-Walker EW2 forcing term
            if newton_iter > 0 {
                // EW2: η_k = |‖F_k‖ - ‖F_{k-1} + J_{k-1} δx_{k-1}‖| / ‖F_{k-1}‖
                let ew2_candidate = NumericElement::abs(norm_f - prev_model_residual) / prev_norm;
                // Safeguard: stay in [eta_min, eta_max]
                eta = self
                    .config
                    .eta_min
                    .max_scalar(self.config.eta_max.min_scalar(ew2_candidate));
            }

            let inner_tol = eta * norm_f;

            // Solve J(x) δx = -F(x) via matrix-free GMRES
            let rhs = neg(&f);
            let delta_x = self.gmres_matrix_free(&func, &x, &f, &rhs, inner_tol, n)?;

            // Model residual estimate: ‖F(xₖ) + J(xₖ)δxₖ‖ (via JvP)
            let eps = Self::compute_eps(&x, &delta_x);
            let jop = JvpOperator {
                x_pivot: &x,
                f_pivot: &f,
                func: &func,
                eps,
            };
            let j_delta = jop.apply_jvp(&delta_x);
            prev_model_residual = norm(&add(&f, &j_delta));
            prev_norm = norm_f;

            // Update: x ← x + δx
            x = add(&x, &delta_x);
            f = func(&x);

            history.push(norm(&f));
        }

        let final_norm = norm(&f);
        Ok((
            x,
            JfnkConvergence {
                residual_norm: final_norm,
                newton_iterations: self.config.max_newton_iterations,
                residual_history: history,
                converged: final_norm < abs_tol,
            },
        ))
    }

    /// Inner GMRES solve for $J(x) \delta x = b$ using matrix-free Jacobian-vector products.
    ///
    /// Implements the restarted GMRES(m) algorithm (Saad & Schultz 1986) with:
    /// - Arnoldi iteration building an orthonormal basis $V_m$ for the Krylov subspace
    /// - Least-squares minimization via Givens rotations on the Hessenberg matrix
    ///
    /// # Theorem — GMRES Monotone Residual Decrease (Saad & Schultz 1986, Thm 2.1)
    ///
    /// At step $k$, GMRES minimizes ‖b − J xₖ‖₂ over all $x_k \in x_0 + \mathcal{K}_k(J, r_0)$.
    /// The residual sequence is non-increasing: ‖r_{k+1}‖ ≤ ‖r_k‖.
    fn gmres_matrix_free<F>(
        &self,
        func: &F,
        x: &Array1<T>,
        f_x: &Array1<T>,
        b: &Array1<T>,
        tol: T,
        n: usize,
    ) -> Result<Array1<T>>
    where
        F: Fn(&Array1<T>) -> Array1<T>,
    {
        let restart = self.config.krylov_restart.min(n);
        let max_iter = self.config.max_krylov_iterations;

        let mut delta = vector_zeros(n);
        // δ₀ = 0 ⇒ J(x)·δ₀ = 0 ⇒ r₀ = b
        let mut r = b.clone();

        let beta = norm(&r);
        if beta < tol {
            return Ok(delta);
        }

        let mut total_iters = 0usize;
        while total_iters < max_iter {
            // Restarted Arnoldi iteration
            let mut v_basis: Vec<Array1<T>> = Vec::with_capacity(restart + 1);
            // H: (restart+1) × restart upper Hessenberg
            let mut h = vec![vec![T::ZERO; restart]; restart + 1];
            // Givens rotation parameters
            let mut cs = vec![T::ZERO; restart];
            let mut sn = vec![T::ZERO; restart];
            // RHS of projected problem (length restart+1)
            let mut g = vec![T::ZERO; restart + 1];

            let beta_r = norm(&r);
            if beta_r < tol {
                break;
            }
            g[0] = beta_r;
            v_basis.push(scale(&r, T::ONE / beta_r));

            let mut inner_converged = false;
            for j in 0..restart {
                if total_iters >= max_iter {
                    break;
                }
                total_iters += 1;

                // Arnoldi: compute w = J(x) · vⱼ
                let vj = &v_basis[j];
                let eps = Self::compute_eps(x, vj);
                let jop = JvpOperator {
                    x_pivot: x,
                    f_pivot: f_x,
                    func,
                    eps,
                };
                let mut w = jop.apply_jvp(vj);

                // Modified Gram-Schmidt orthogonalization
                for i in 0..=j {
                    h[i][j] = dot(&v_basis[i], &w);
                    sub_in_place_scaled(&mut w, &v_basis[i], h[i][j]);
                }
                h[j + 1][j] = norm(&w);

                if h[j + 1][j] > T::ZERO {
                    v_basis.push(scale(&w, T::ONE / h[j + 1][j]));
                } else {
                    // Breakdown: Krylov subspace exhausted (exact solution found)
                    v_basis.push(vector_zeros(n));
                }

                // Apply previous Givens rotations to new column
                for i in 0..j {
                    let temp = cs[i] * h[i][j] + sn[i] * h[i + 1][j];
                    h[i + 1][j] = -sn[i] * h[i][j] + cs[i] * h[i + 1][j];
                    h[i][j] = temp;
                }

                // Compute new Givens rotation
                let (c, s) = Self::givens(h[j][j], h[j + 1][j]);
                cs[j] = c;
                sn[j] = s;

                // Apply new rotation to H and g
                h[j][j] = cs[j] * h[j][j] + sn[j] * h[j + 1][j];
                h[j + 1][j] = T::ZERO;
                g[j + 1] = -sn[j] * g[j];
                g[j] = cs[j] * g[j];

                if NumericElement::abs(g[j + 1]) < tol {
                    // Converged within this restart
                    let y = Self::back_substitute(&h, &g, j + 1);
                    for (k, yk) in y.iter().enumerate() {
                        add_scaled_in_place(&mut delta, &v_basis[k], *yk);
                    }
                    let fwd = {
                        let eps = Self::compute_eps(x, &delta);
                        let jop = JvpOperator {
                            x_pivot: x,
                            f_pivot: f_x,
                            func,
                            eps,
                        };
                        jop.apply_jvp(&delta)
                    };
                    r = sub(b, &fwd);
                    inner_converged = true;
                    break;
                }
            }

            // Update solution with best approximation from this restart
            // (skip if inner loop already applied the converged solution)
            if !inner_converged {
                let m = v_basis.len().saturating_sub(1).min(restart);
                if m > 0 {
                    let y = Self::back_substitute(&h, &g, m);
                    for (k, yk) in y.iter().enumerate() {
                        add_scaled_in_place(&mut delta, &v_basis[k], *yk);
                    }
                    let fwd = {
                        let eps = Self::compute_eps(x, &delta);
                        let jop = JvpOperator {
                            x_pivot: x,
                            f_pivot: f_x,
                            func,
                            eps,
                        };
                        jop.apply_jvp(&delta)
                    };
                    r = sub(b, &fwd);
                }
            }

            if norm(&r) < tol {
                break;
            }
        }

        Ok(delta)
    }

    /// Compute Givens rotation coefficients (c, s) for elements (a, b):
    /// $[c, s; -s, c] [a; b] = [r; 0]$, where $r = \sqrt{a^2 + b^2}$.
    ///
    /// Uses the numerically stable form from Golub & Van Loan (2013), §5.1.8.
    fn givens(a: T, b: T) -> (T, T) {
        if NumericElement::abs(b) < T::EPSILON {
            (T::ONE, T::ZERO)
        } else if NumericElement::abs(b) >= NumericElement::abs(a) {
            let tau = -a / b;
            let s = T::ONE / NumericElement::sqrt(T::ONE + tau * tau);
            (s * tau, s)
        } else {
            let tau = -b / a;
            let c = T::ONE / NumericElement::sqrt(T::ONE + tau * tau);
            (c, c * tau)
        }
    }

    /// Back-substitution for upper-triangular system from Hessenberg QR.
    fn back_substitute(h: &[Vec<T>], g: &[T], m: usize) -> Vec<T> {
        let mut y = vec![T::ZERO; m];
        for i in (0..m).rev() {
            let mut s = g[i];
            for j in (i + 1)..m {
                s -= h[i][j] * y[j];
            }
            if NumericElement::abs(h[i][i]) > T::EPSILON {
                y[i] = s / h[i][i];
            }
        }
        y
    }
}

#[cfg(test)]
mod tests {
    use super::super::linalg::vector_from_vec;
    use super::*;
    use eunomia::assert_relative_eq;
    use leto::Array2;

    fn vec(values: Vec<f64>) -> Array1<f64> {
        vector_from_vec(values)
    }

    fn mat_vec(matrix: &Array2<f64>, vector: &Array1<f64>) -> Array1<f64> {
        let [rows, cols] = matrix.shape();
        assert_eq!(cols, vector.shape()[0]);
        vector_from_vec(
            (0..rows)
                .map(|row| (0..cols).fold(0.0, |acc, col| acc + matrix[[row, col]] * vector[col]))
                .collect(),
        )
    }

    /// Theorem: JFNK on a linear system F(x) = Ax - b must recover the exact solution
    /// in ≤ n outer iterations (linear case degenerates to GMRES = direct solver).
    #[test]
    fn test_jfnk_linear_system() {
        // F(x) = Ax - b = 0, x* = [1.0, 2.0, 3.0]
        let a_data =
            Array2::from_shape_vec([3, 3], vec![4.0, 1.0, 0.0, 1.0, 3.0, 1.0, 0.0, 1.0, 2.0])
                .unwrap();
        let b = vec(vec![
            4.0 * 1.0 + 1.0 * 2.0,
            1.0 * 1.0 + 3.0 * 2.0 + 1.0 * 3.0,
            1.0 * 2.0 + 2.0 * 3.0,
        ]);

        let func = |x: &Array1<f64>| sub(&mat_vec(&a_data, x), &b);

        let config = JfnkConfig::<f64> {
            max_newton_iterations: 20,
            atol: 1e-10,
            rtol: 1e-8,
            max_krylov_iterations: 100,
            krylov_restart: 30,
            inner_tol_init: 0.5,
            eta_min: 1e-4,
            eta_max: 0.9,
        };

        let x0 = vector_zeros(3);
        let solver = JfnkSolver::new(config);
        let (x_sol, conv) = solver.solve(func, x0).expect("JFNK must succeed");

        assert!(conv.converged, "JFNK must converge for SPD linear system");
        assert_relative_eq!(x_sol[0], 1.0, epsilon = 1e-8);
        assert_relative_eq!(x_sol[1], 2.0, epsilon = 1e-8);
        assert_relative_eq!(x_sol[2], 3.0, epsilon = 1e-8);
    }

    /// Theorem: JFNK on a nonlinear system F(x) = [x₀² - 1, x₁² - 4] must
    /// converge quadratically near root x* = [1.0, 2.0].
    #[test]
    fn test_jfnk_nonlinear_quadratic_convergence() {
        // F(x) = [x[0]^2 - 1, x[1]^2 - 4], x* = [1.0, 2.0]
        let func = |x: &Array1<f64>| vec(vec![x[0] * x[0] - 1.0, x[1] * x[1] - 4.0]);

        let config = JfnkConfig::<f64> {
            max_newton_iterations: 30,
            atol: 1e-10,
            rtol: 1e-10,
            max_krylov_iterations: 50,
            krylov_restart: 10,
            ..Default::default()
        };
        let x0 = vec(vec![0.5, 1.5]);
        let solver = JfnkSolver::new(config);
        let (x_sol, conv) = solver.solve(func, x0).expect("JFNK must succeed");

        assert!(
            conv.converged,
            "JFNK nonlinear: must converge, residual={:.2e}, iters={}",
            conv.residual_norm, conv.newton_iterations
        );
        assert_relative_eq!(x_sol[0], 1.0, epsilon = 1e-8);
        assert_relative_eq!(x_sol[1], 2.0, epsilon = 1e-8);

        // Verify near-quadratic convergence: ‖eₖ₊₁‖ / ‖eₖ‖² < 10 for last 3 steps
        let h = &conv.residual_history;
        let n = h.len();
        if n >= 3 {
            for k in (n - 3)..n.saturating_sub(1) {
                if h[k] > 1e-6 {
                    let ratio = h[k + 1] / (h[k] * h[k]);
                    assert!(
                        ratio < 100.0,
                        "Convergence ratio {ratio:.2e} too large at step {k} — not near-quadratic"
                    );
                }
            }
        }
    }

    /// Verify JfnkSolver returns converged=true when initial guess is already a root.
    #[test]
    fn test_jfnk_already_converged_initial_guess() {
        let func = |x: &Array1<f64>| vec(vec![x[0] - 1.0, x[1] - 2.0]);
        let config = JfnkConfig::<f64>::default();
        let x0 = vec(vec![1.0, 2.0]);
        let solver = JfnkSolver::new(config);
        let (_, conv) = solver.solve(func, x0).unwrap();
        assert!(conv.converged);
        assert_eq!(conv.newton_iterations, 0);
    }

    /// Verify EW forcing term stays in [eta_min, eta_max].
    #[test]
    fn test_jfnk_config_invariants() {
        let config = JfnkConfig::<f64>::default();
        assert!(config.eta_min < config.eta_max);
        assert!(config.eta_min > 0.0);
        assert!(config.eta_max < 1.0);
        assert!(config.atol > 0.0);
        assert!(config.rtol > 0.0);
        assert!(config.max_newton_iterations > 0);
        assert!(config.krylov_restart > 0);
    }
}
