//! Anderson Acceleration (Extrapolation) for Fixed-Point Iterations
//!
//! # Theorem — Anderson Acceleration Convergence (Anderson, 1965)
//!
//! Let $G: \mathbb{R}^n \to \mathbb{R}^n$ be a contractive mapping with fixed point
//! $x^* = G(x^*)$. Anderson Acceleration computes an accelerated update:
//!
//! $$ \mathbf{x}_{k+1} = \mathbf{x}_k + \beta \mathbf{f}_k - (\Delta \mathbf{X}_k + \beta \Delta \mathbf{F}_k) \gamma $$
//!
//! where $\mathbf{f}_k = G(\mathbf{x}_k) - \mathbf{x}_k$ is the residual, $\Delta \mathbf{X}_k$ and
//! $\Delta \mathbf{F}_k$ are matrices of the last $m$ step differences, and $\gamma$ solves:
//!
//! $$ \min_\gamma \|\mathbf{f}_k - \Delta \mathbf{F}_k \gamma\|_2 $$
//!
//! Locally achieves superlinear convergence without an explicit Jacobian.
//!
//! # Theorem — MGS-QR Anderson vs Normal Equations (Walker & Ni 2011, Thm 2.1)
//!
//! When solving $\min_\gamma \|f - \Delta F \gamma\|_2$, two approaches are possible:
//!
//! **Normal equations** (Type-I): $({\Delta F}^T \Delta F)\gamma = {\Delta F}^T f$.
//! - Condition number: $\kappa(\Delta F^T \Delta F) = \kappa(\Delta F)^2$
//! - Numerically unstable when $\Delta F$ columns are nearly linearly dependent.
//!
//! **QR factorization** (Type-II): $\Delta F = QR \Rightarrow R\gamma = Q^T f$.
//! - Condition number: $\kappa(R) = \kappa(\Delta F)$
//! - Stable: MGS-QR halves the sensitivity to near-linear-dependence in history.
//!
//! **Proof sketch**: For $\Delta F = QR$ (thin QR), the unique least-squares solution
//! is $\gamma^* = R^{-1} Q^T f$ whenever $\Delta F$ has full column rank. The
//! backward-stable MGS process produces $\|Q^T Q - I\| = O(\epsilon_{\rm mach} \kappa(\Delta F))$.
//! The normal equations approach amplifies this error by $\kappa(\Delta F)$, giving
//! $O(\epsilon_{\rm mach} \kappa(\Delta F)^2)$ rounding error in $\gamma^*$.
//!
//! **Reference**: Walker, H.F. & Ni, P. (2011). Anderson acceleration for fixed-point
//! iterations. *SIAM J. Numer. Anal.* 49(4):1715–1735.
//!
//! # Theorem — VecDeque O(1) history eviction (GAP-PERF-004)
//!
//! History eviction via `Vec::remove(0)` is O(m) (memory shift of m vectors).
//! `VecDeque::pop_front()` is O(1) (pointer rotation on ring buffer).
//! For history depth m=5 and 10³ outer iterations, total shift cost drops from
//! O(5 × 10³) = 5000 ops to O(10³) = 1000 ops in pointer increments.

use nalgebra::{DMatrix, DVector, RealField};
use num_traits::{Float, FromPrimitive};
use std::collections::VecDeque;

/// Method for solving the Anderson least-squares subproblem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AndersonMethod {
    /// Solve via normal equations (ΔFᵀΔF)γ = ΔFᵀf using SVD.
    /// Simple but condition number = κ(ΔF)².
    NormalEquations,
    /// Solve via incremental Modified Gram-Schmidt QR factorization.
    /// Condition number = κ(ΔF). Recommended for ill-conditioned histories.
    #[default]
    QR,
}

/// Configuration for Anderson Acceleration
#[derive(Debug, Clone)]
pub struct AndersonConfig<T: RealField> {
    /// Maximum history size ($m$)
    pub history_depth: usize,
    /// Relaxation / damping factor ($\beta \in (0, 1]$)
    pub relaxation: T,
    /// Minimum tolerance for SVD or QR to prevent singular matrix issues in least-squares solve
    pub drop_tolerance: T,
    /// Algorithm used to solve the Anderson least-squares subproblem.
    pub method: AndersonMethod,
}

impl<T: RealField + FromPrimitive> Default for AndersonConfig<T> {
    fn default() -> Self {
        Self {
            history_depth: 5,
            relaxation: T::from_f64(1.0).unwrap_or_else(T::one),
            drop_tolerance: T::from_f64(1e-12).unwrap_or_else(T::zero),
            method: AndersonMethod::QR,
        }
    }
}

/// Internal QR state for incremental Modified Gram-Schmidt updates.
///
/// Maintains a thin QR factorization ΔF = Q R where:
/// - `q_cols`: columns of Q (orthonormal basis for col-span of ΔF), stored as Vec<DVector>
/// - `r_mat`:  upper-triangular R, stored as a dense DMatrix (m×m, m small)
///
/// Columns are appended when new history arrives and dropped (oldest first)
/// to maintain `history_depth` limit — using a simple shift on the small (m≤10) R.
#[derive(Debug)]
struct QrState<T: RealField + Copy> {
    /// Orthonormal columns of Q (N-dimensional), bounded to `m` entries.
    q_cols: VecDeque<DVector<T>>,
    /// Upper-triangular R stored as dense column-major m×m.
    /// `r_mat[(i, j)]` = R[i,j], j-th column of R corresponds to j-th ΔF column.
    r_mat: DMatrix<T>,
    /// Maximum history depth m.
    max_depth: usize,
    /// Drop tolerance for near-zero diagonal entries in R.
    drop_tol: T,
}

impl<T: RealField + Copy + Float + std::fmt::Debug> QrState<T> {
    fn new(max_depth: usize, drop_tol: T) -> Self {
        Self {
            q_cols: VecDeque::with_capacity(max_depth + 1),
            r_mat: DMatrix::zeros(0, 0),
            max_depth,
            drop_tol,
        }
    }

    /// Append a new column `v` to ΔF, extending Q and R by one column via MGS.
    ///
    /// If history is full, drop the oldest column first (O(m²) R shift).
    fn append_column(&mut self, v: &DVector<T>, _m_hint: usize) {
        let m_current = self.q_cols.len();
        // Drop oldest column when at capacity
        if m_current >= self.max_depth && !self.q_cols.is_empty() {
            // Remove leftmost q column (O(m) ring-buffer pop)
            self.q_cols.pop_front();
            // Shift R: remove first column and first row (O(m²) small matrix op)
            let m = self.q_cols.len();
            if m > 0 {
                // R shrinks from (m+1)×(m+1) to m×m by dropping row 0 and col 0
                let new_r = self.r_mat.view((1, 1), (m, m)).into_owned();
                self.r_mat = new_r;
            } else {
                self.r_mat = DMatrix::zeros(0, 0);
            }
        }

        let m_new = self.q_cols.len(); // columns already in Q after potential drop

        // MGS orthogonalization: project `v` against all existing Q columns
        let mut w = v.clone();
        let mut r_col = DVector::zeros(m_new + 1);
        for (j, qj) in self.q_cols.iter().enumerate() {
            let alpha = qj.dot(&w);
            r_col[j] = alpha;
            w -= qj * alpha;
        }

        // Norm of residual = R[m_new, m_new]
        let norm_w = Float::sqrt(w.dot(&w));
        r_col[m_new] = norm_w;

        // Append new Q column (or handle near-zero case by skipping)
        if norm_w > self.drop_tol {
            self.q_cols.push_back(w / norm_w);

            // Expand R by one column on the right and one row at the bottom
            let new_m = m_new + 1;
            let mut new_r = DMatrix::zeros(new_m, new_m);
            if m_new > 0 && self.r_mat.nrows() == m_new && self.r_mat.ncols() == m_new {
                new_r
                    .view_mut((0, 0), (m_new, m_new))
                    .copy_from(&self.r_mat);
            }
            for i in 0..=m_new {
                new_r[(i, m_new)] = r_col[i];
            }
            self.r_mat = new_r;
        }
        // If norm_w ≤ drop_tol: column is linearly dependent; discard silently.
    }

    /// Solve the least-squares problem min ‖f − ΔF γ‖₂ = min ‖f − QR γ‖₂
    /// via back-substitution: γ = R⁻¹ (Qᵀ f).
    ///
    /// Returns `None` if R is degenerate (any diagonal < drop_tol).
    fn solve_least_squares(&self, f: &DVector<T>) -> Option<DVector<T>> {
        let m = self.q_cols.len();
        if m == 0 || self.r_mat.nrows() != m || self.r_mat.ncols() != m {
            return None;
        }

        // Compute rhs = Qᵀ f  (m-vector)
        let mut rhs = DVector::zeros(m);
        for (j, qj) in self.q_cols.iter().enumerate() {
            rhs[j] = qj.dot(f);
        }

        // Back-substitution: R γ = rhs (R is upper-triangular m×m)
        let mut gamma = DVector::zeros(m);
        for i in (0..m).rev() {
            let diag = self.r_mat[(i, i)];
            if Float::abs(diag) < self.drop_tol {
                return None; // Degenerate column — history will be cleared by caller.
            }
            let mut s = rhs[i];
            for j in (i + 1)..m {
                s -= self.r_mat[(i, j)] * gamma[j];
            }
            gamma[i] = s / diag;
        }
        Some(gamma)
    }

    /// Clear all QR state.
    fn reset(&mut self) {
        self.q_cols.clear();
        self.r_mat = DMatrix::zeros(0, 0);
    }
}

/// Anderson Accelerator for non-linear Picard / Fixed-Point iterations.
///
/// Supports two subproblem methods via `AndersonConfig::method`:
/// - `QR` (default): incremental MGS-QR, condition number κ(ΔF)
/// - `NormalEquations`: SVD on normal equations, condition number κ(ΔF)²
#[derive(Debug)]
pub struct AndersonAccelerator<T: RealField + Copy> {
    config: AndersonConfig<T>,

    /// Previous state vector $x_{k-1}$
    prev_x: Option<DVector<T>>,
    /// Previous residual vector $f_{k-1} = G(x_{k-1}) - x_{k-1}$
    prev_f: Option<DVector<T>>,

    /// History of state differences: $\Delta X = [ \Delta x_{k-m}, \dots, \Delta x_{k-1} ]$
    /// `VecDeque` for O(1) front eviction (GAP-PERF-004).
    delta_x: VecDeque<DVector<T>>,
    /// History of residual differences: $\Delta F = [ \Delta f_{k-m}, \dots, \Delta f_{k-1} ]$
    /// `VecDeque` for O(1) front eviction (GAP-PERF-004).
    delta_f: VecDeque<DVector<T>>,

    /// Incremental QR state (used only when method == QR).
    qr_state: Option<QrState<T>>,
}

impl<T: RealField + Copy + Float + std::fmt::Debug> AndersonAccelerator<T> {
    /// Create a new Anderson Acceleration context.
    #[must_use]
    pub fn new(config: AndersonConfig<T>) -> Self {
        let depth = config.history_depth;
        let qr_state = if config.method == AndersonMethod::QR {
            Some(QrState::new(config.history_depth, config.drop_tolerance))
        } else {
            None
        };
        Self {
            config,
            prev_x: None,
            prev_f: None,
            delta_x: VecDeque::with_capacity(depth + 1),
            delta_f: VecDeque::with_capacity(depth + 1),
            qr_state,
        }
    }

    /// Calculate the next accelerated state $\mathbf{x}_{k+1}$.
    ///
    /// # Arguments
    /// * `x`   — Current input state vector $\mathbf{x}_k$.
    /// * `g_x` — Output of fixed-point operator $G(\mathbf{x}_k)$.
    ///
    /// # Returns
    /// Accelerated next state $\mathbf{x}_{k+1}$.
    #[allow(clippy::many_single_char_names)]
    pub fn compute_next(&mut self, x: &DVector<T>, g_x: &DVector<T>) -> DVector<T> {
        let f = g_x - x;

        let mut x_next = x + &f * self.config.relaxation;

        if let (Some(prev_x), Some(prev_f)) = (&self.prev_x, &self.prev_f) {
            let dx = x - prev_x;
            let df = &f - prev_f;

            // Evict oldest history entry if at capacity — O(1) for VecDeque (GAP-PERF-004)
            if self.delta_x.len() >= self.config.history_depth {
                self.delta_x.pop_front();
                self.delta_f.pop_front();
            }
            let m_before = self.delta_f.len();
            self.delta_x.push_back(dx);
            self.delta_f.push_back(df.clone());

            // Solve Anderson subproblem
            let gamma_opt = match self.config.method {
                AndersonMethod::QR => {
                    // Update QR factorization with new ΔF column, then solve
                    if let Some(qr) = &mut self.qr_state {
                        qr.append_column(&df, m_before);
                        qr.solve_least_squares(&f)
                    } else {
                        None
                    }
                }
                AndersonMethod::NormalEquations => self.solve_normal_equations(&f),
            };

            match gamma_opt {
                Some(gamma) => {
                    // x_next = x + β·f − (ΔX + β·ΔF) · γ
                    for (j, g_j) in gamma.iter().enumerate() {
                        let term = &self.delta_x[j] + &self.delta_f[j] * self.config.relaxation;
                        x_next -= term * *g_j;
                    }
                }
                None => {
                    // Degenerate subproblem — clear history and recover with plain relaxation
                    self.delta_x.clear();
                    self.delta_f.clear();
                    if let Some(qr) = &mut self.qr_state {
                        qr.reset();
                    }
                }
            }
        }

        self.prev_x = Some(x.clone());
        self.prev_f = Some(f);

        x_next
    }

    /// Solve via normal equations using SVD (original Type-I method).
    fn solve_normal_equations(&self, f: &DVector<T>) -> Option<DVector<T>> {
        let m = self.delta_f.len();
        if m == 0 {
            return None;
        }
        let n = f.len();
        let mut mat_df = DMatrix::zeros(n, m);
        for (j, df_col) in self.delta_f.iter().enumerate() {
            mat_df.set_column(j, df_col);
        }
        let df_t_df = mat_df.transpose() * &mat_df;
        let df_t_f = mat_df.transpose() * f;
        let svd = df_t_df.svd(true, true);
        svd.solve(&df_t_f, self.config.drop_tolerance).ok()
    }

    /// Reset the acceleration history (useful if the solver detects stagnation).
    pub fn reset(&mut self) {
        self.prev_x = None;
        self.prev_f = None;
        self.delta_x.clear();
        self.delta_f.clear();
        if let Some(qr) = &mut self.qr_state {
            qr.reset();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::DVector;

    fn make_linear_system() -> (nalgebra::DMatrix<f64>, DVector<f64>) {
        // G(x) = A x + b, fixed point at x* = [1.0, 2.0]
        let a = nalgebra::DMatrix::from_row_slice(2, 2, &[0.5, 0.1, 0.0, 0.5]);
        let b = DVector::from_vec(vec![1.0 - (0.5 * 1.0 + 0.1 * 2.0), 2.0 - 0.5 * 2.0]);
        (a, b)
    }

    fn run_fixed_point(config: AndersonConfig<f64>, iters: usize) -> DVector<f64> {
        let (a, b) = make_linear_system();
        let mut accelerator = AndersonAccelerator::new(config);
        let mut x = DVector::from_vec(vec![0.0, 0.0]);
        for _ in 0..iters {
            let g_x = &a * &x + &b;
            x = accelerator.compute_next(&x, &g_x);
        }
        x
    }

    /// Verify QR variant converges to x* = [1.0, 2.0] within 1e-8.
    #[test]
    fn test_anderson_qr_convergence() {
        let config = AndersonConfig::<f64> {
            history_depth: 3,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::QR,
        };
        let x = run_fixed_point(config, 30);
        assert_relative_eq!(x[0], 1.0, epsilon = 1e-8);
        assert_relative_eq!(x[1], 2.0, epsilon = 1e-8);
    }

    /// QR and NormalEquations variants must converge to the same fixed point.
    #[test]
    fn test_anderson_qr_equivalence_to_normal_equations() {
        let config_qr = AndersonConfig::<f64> {
            history_depth: 3,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::QR,
        };
        let config_ne = AndersonConfig::<f64> {
            method: AndersonMethod::NormalEquations,
            ..config_qr.clone()
        };
        let x_qr = run_fixed_point(config_qr, 30);
        let x_ne = run_fixed_point(config_ne, 30);
        assert_relative_eq!(x_qr[0], x_ne[0], epsilon = 1e-7);
        assert_relative_eq!(x_qr[1], x_ne[1], epsilon = 1e-7);
    }

    /// VecDeque history must not exceed m entries after m+5 steps (GAP-PERF-004).
    #[test]
    fn test_vecdeque_history_bounded() {
        let m = 3usize;
        let config = AndersonConfig::<f64> {
            history_depth: m,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::QR,
        };
        let (a, b) = make_linear_system();
        let mut acc = AndersonAccelerator::new(config);
        let mut x = DVector::from_vec(vec![0.0, 0.0]);
        for _ in 0..(m + 5) {
            let g_x = &a * &x + &b;
            x = acc.compute_next(&x, &g_x);
            assert!(
                acc.delta_x.len() <= m,
                "delta_x history len {} exceeds m={}",
                acc.delta_x.len(),
                m
            );
            assert!(
                acc.delta_f.len() <= m,
                "delta_f history len {} exceeds m={}",
                acc.delta_f.len(),
                m
            );
        }
    }

    /// Original NormalEquations convergence test preserved (regression guard).
    #[test]
    fn test_anderson_acceleration_linear_convergence() {
        let config = AndersonConfig::<f64> {
            history_depth: 3,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::NormalEquations,
        };
        let x = run_fixed_point(config, 20);
        assert!((x[0] - 1.0).abs() < 1e-6);
        assert!((x[1] - 2.0).abs() < 1e-6);
    }

    /// QR should handle near-collinear history without divergence (ill-conditioning test).
    /// When ΔF columns are nearly identical, the normal equations approach is prone to
    /// numerical rank deficiency — QR should silently discard bad columns.
    #[test]
    fn test_anderson_qr_ill_conditioned_history_graceful() {
        // Very slowly-converging fixed point: contraction ratio 0.99 → lots of similar steps
        let a = nalgebra::DMatrix::from_row_slice(2, 2, &[0.99, 0.0, 0.0, 0.99]);
        let b = DVector::from_vec(vec![0.01, 0.02]); // x* = [1.0, 2.0]

        let config = AndersonConfig::<f64> {
            history_depth: 5,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::QR,
        };
        let mut acc = AndersonAccelerator::new(config);
        let mut x = DVector::from_vec(vec![0.0, 0.0]);

        // Run until convergence or 200 iterations — must not panic or diverge
        for _ in 0..200 {
            let g_x = &a * &x + &b;
            x = acc.compute_next(&x, &g_x);
            if (x[0] - 1.0).abs() < 1e-6 && (x[1] - 2.0).abs() < 1e-6 {
                break;
            }
        }
        assert!(
            (x[0] - 1.0).abs() < 1e-4,
            "QR Anderson must converge on ill-conditioned problem: x[0]={:.6}",
            x[0]
        );
        assert!(
            (x[1] - 2.0).abs() < 1e-4,
            "QR Anderson must converge on ill-conditioned problem: x[1]={:.6}",
            x[1]
        );
    }
}
