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

use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};
use std::collections::VecDeque;

use super::linalg::{
    add_scaled, add_scaled_in_place, dot, matrix_zeros, scale, sub, vector_len, vector_zeros,
};

/// Method for solving the Anderson least-squares subproblem.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AndersonMethod {
    /// Solve via normal equations (ΔFᵀΔF)γ = ΔFᵀf using a pivoted dense solve.
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
    /// Minimum tolerance for QR or pivoted solves to prevent singular matrix issues
    /// in least-squares solve.
    pub drop_tolerance: T,
    /// Algorithm used to solve the Anderson least-squares subproblem.
    pub method: AndersonMethod,
}

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

impl<T: RealField + FloatElement> Default for AndersonConfig<T> {
    fn default() -> Self {
        Self {
            history_depth: 5,
            relaxation: from_f64(1.0),
            drop_tolerance: from_f64(1e-12),
            method: AndersonMethod::QR,
        }
    }
}

/// Internal QR state for incremental Modified Gram-Schmidt updates.
///
/// Maintains a thin QR factorization ΔF = Q R where:
/// - `q_cols`: columns of Q (orthonormal basis for col-span of ΔF), stored as Leto `Array1`
/// - `r_mat`:  upper-triangular R, stored as a dense Leto `Array2` (m×m, m small)
///
/// Columns are appended when new history arrives and dropped (oldest first)
/// to maintain `history_depth` limit — using a simple shift on the small (m≤10) R.
#[derive(Debug)]
struct QrState<T: RealField + Copy> {
    /// Orthonormal columns of Q (N-dimensional), bounded to `m` entries.
    q_cols: VecDeque<Array1<T>>,
    /// Upper-triangular R stored as dense column-major m×m.
    /// `r_mat[(i, j)]` = R[i,j], j-th column of R corresponds to j-th ΔF column.
    r_mat: Array2<T>,
    /// Maximum history depth m.
    max_depth: usize,
    /// Drop tolerance for near-zero diagonal entries in R.
    drop_tol: T,
}

impl<T: RealField + Copy + FloatElement + std::fmt::Debug> QrState<T> {
    fn new(max_depth: usize, drop_tol: T) -> Self {
        Self {
            q_cols: VecDeque::with_capacity(max_depth + 1),
            r_mat: matrix_zeros(0, 0),
            max_depth,
            drop_tol,
        }
    }

    /// Append a new column `v` to ΔF, extending Q and R by one column via MGS.
    ///
    /// If the history is full, the oldest column is dropped first (O(m²) R
    /// shift). When the new direction is linearly dependent on the active
    /// history (norm_w ≤ `drop_tol`), the column is **silently rejected** and
    /// the QR state remains unchanged.
    ///
    /// # Returns
    ///
    /// `true` iff the column was accepted (Q/R grown by one); `false` iff it
    /// was rejected as near-collinear. Callers that keep a parallel
    /// `delta_x`/`delta_f` history **must** match the result: on `false`, the
    /// just-pushed `(Δx, Δf)` pair must be evicted from the parallel history
    /// to preserve the lockstep invariant
    /// `self.q_cols.len() == delta_x.len() == delta_f.len()`.
    fn append_column(&mut self, v: &Array1<T>, _m_hint: usize) -> bool {
        let m_current = self.q_cols.len();
        // Drop oldest column when at capacity
        if m_current >= self.max_depth && !self.q_cols.is_empty() {
            // Remove leftmost q column (O(m) ring-buffer pop)
            self.q_cols.pop_front();
            // Shift R: remove first column and first row (O(m²) small matrix op)
            let m = self.q_cols.len();
            if m > 0 {
                // R shrinks from (m+1)×(m+1) to m×m by dropping row 0 and col 0
                let mut new_r = matrix_zeros(m, m);
                for row in 0..m {
                    for col in 0..m {
                        new_r[[row, col]] = self.r_mat[[row + 1, col + 1]];
                    }
                }
                self.r_mat = new_r;
            } else {
                self.r_mat = matrix_zeros(0, 0);
            }
        }

        let m_new = self.q_cols.len(); // columns already in Q after potential drop

        // MGS orthogonalization: project `v` against all existing Q columns
        let mut w = v.clone();
        let mut r_col = vector_zeros(m_new + 1);
        for (j, qj) in self.q_cols.iter().enumerate() {
            let alpha = dot(qj, &w);
            r_col[[j]] = alpha;
            add_scaled_in_place(&mut w, qj, T::ZERO - alpha);
        }

        // Norm of residual = R[m_new, m_new]
        let norm_w = NumericElement::sqrt(dot(&w, &w));
        r_col[[m_new]] = norm_w;

        // Append new Q column, or signal rejection if linearly dependent.
        if norm_w <= self.drop_tol {
            // Column is linearly dependent on the active history — discard.
            // Returning `false` lets the caller evict the corresponding
            // entry from the parallel `delta_x`/`delta_f` deques so the
            // lockstep invariant holds.
            return false;
        }

        self.q_cols.push_back(scale(&w, T::ONE / norm_w));

        // Expand R by one column on the right and one row at the bottom
        let new_m = m_new + 1;
        let mut new_r = matrix_zeros(new_m, new_m);
        if m_new > 0 && self.r_mat.shape() == [m_new, m_new] {
            for row in 0..m_new {
                for col in 0..m_new {
                    new_r[[row, col]] = self.r_mat[[row, col]];
                }
            }
        }
        for i in 0..=m_new {
            new_r[[i, m_new]] = r_col[[i]];
        }
        self.r_mat = new_r;

        true
    }

    /// Solve the least-squares problem min ‖f − ΔF γ‖₂ = min ‖f − QR γ‖₂
    /// via back-substitution: γ = R⁻¹ (Qᵀ f).
    ///
    /// Returns `None` if R is degenerate (any diagonal < drop_tol).
    fn solve_least_squares(&self, f: &Array1<T>) -> Option<Array1<T>> {
        let m = self.q_cols.len();
        if m == 0 || self.r_mat.shape() != [m, m] {
            return None;
        }

        // Compute rhs = Qᵀ f  (m-vector)
        let mut rhs = vector_zeros(m);
        for (j, qj) in self.q_cols.iter().enumerate() {
            rhs[[j]] = dot(qj, f);
        }

        // Back-substitution: R γ = rhs (R is upper-triangular m×m)
        let mut gamma = vector_zeros(m);
        for i in (0..m).rev() {
            let diag = self.r_mat[[i, i]];
            if NumericElement::abs(diag) < self.drop_tol {
                return None; // Degenerate column — history will be cleared by caller.
            }
            let mut s = rhs[[i]];
            for j in (i + 1)..m {
                s -= self.r_mat[[i, j]] * gamma[[j]];
            }
            gamma[[i]] = s / diag;
        }
        Some(gamma)
    }

    /// Clear all QR state.
    fn reset(&mut self) {
        self.q_cols.clear();
        self.r_mat = matrix_zeros(0, 0);
    }
}

/// Anderson Accelerator for non-linear Picard / Fixed-Point iterations.
///
/// Supports two subproblem methods via `AndersonConfig::method`:
/// - `QR` (default): incremental MGS-QR, condition number κ(ΔF)
/// - `NormalEquations`: pivoted normal equations solve, condition number κ(ΔF)²
#[derive(Debug)]
pub struct AndersonAccelerator<T: RealField + Copy> {
    config: AndersonConfig<T>,

    /// Previous state vector $x_{k-1}$
    prev_x: Option<Array1<T>>,
    /// Previous residual vector $f_{k-1} = G(x_{k-1}) - x_{k-1}$
    prev_f: Option<Array1<T>>,

    /// History of state differences: $\Delta X = [ \Delta x_{k-m}, \dots, \Delta x_{k-1} ]$
    /// `VecDeque` for O(1) front eviction (GAP-PERF-004).
    delta_x: VecDeque<Array1<T>>,
    /// History of residual differences: $\Delta F = [ \Delta f_{k-m}, \dots, \Delta f_{k-1} ]$
    /// `VecDeque` for O(1) front eviction (GAP-PERF-004).
    delta_f: VecDeque<Array1<T>>,

    /// Incremental QR state (used only when method == QR).
    qr_state: Option<QrState<T>>,
}

impl<T: RealField + Copy + FloatElement + std::fmt::Debug> AndersonAccelerator<T> {
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
    pub fn compute_next(&mut self, x: &Array1<T>, g_x: &Array1<T>) -> Array1<T> {
        let f = sub(g_x, x);

        let mut x_next = add_scaled(x, &f, self.config.relaxation);

        if let (Some(prev_x), Some(prev_f)) = (&self.prev_x, &self.prev_f) {
            let dx = sub(x, prev_x);
            let df = sub(&f, prev_f);

            // Evict oldest history entry if at capacity — O(1) for VecDeque (GAP-PERF-004)
            if self.delta_x.len() >= self.config.history_depth {
                self.delta_x.pop_front();
                self.delta_f.pop_front();
            }
            self.delta_x.push_back(dx);
            self.delta_f.push_back(df.clone());

            // Solve Anderson subproblem.
            //
            // QR-path lockstep invariant:
            //     delta_x.len() == delta_f.len() == qr.q_cols.len()
            // holds at every iteration. When the QR silently rejects a
            // near-collinear column (norm_w ≤ drop_tol), the just-pushed
            // pair `(Δx, Δf)` is evicted from the deques so future
            // `compute_next` calls index an in-sync (γ, Δx, Δf) tuple.
            let gamma_opt = match self.config.method {
                AndersonMethod::QR => {
                    let accepted = if let Some(qr) = &mut self.qr_state {
                        qr.append_column(&df, 0)
                    } else {
                        false
                    };
                    if accepted {
                        if let Some(qr) = &self.qr_state {
                            qr.solve_least_squares(&f)
                        } else {
                            None
                        }
                    } else {
                        // QR rejected this column — drop the parallel
                        // deque entry to keep the indices aligned.
                        self.delta_x.pop_back();
                        self.delta_f.pop_back();
                        None
                    }
                }
                AndersonMethod::NormalEquations => self.solve_normal_equations(&f),
            };

            // Invariant (debug build only): the QR path must keep
            // delta_x, delta_f, and qr.q_cols in lockstep.
            #[cfg(debug_assertions)]
            if let Some(qr) = &self.qr_state {
                debug_assert_eq!(
                    self.delta_x.len(),
                    qr.q_cols.len(),
                    "Anderson QR lockstep violated: \
                     delta_x.len()={}, qr.q_cols.len()={}",
                    self.delta_x.len(),
                    qr.q_cols.len(),
                );
                debug_assert_eq!(
                    self.delta_f.len(),
                    qr.q_cols.len(),
                    "Anderson QR lockstep violated: \
                     delta_f.len()={}, qr.q_cols.len()={}",
                    self.delta_f.len(),
                    qr.q_cols.len(),
                );
            }

            if let Some(gamma) = gamma_opt {
                // x_next = x + β·f − (ΔX + β·ΔF) · γ
                for j in 0..vector_len(&gamma) {
                    let term =
                        add_scaled(&self.delta_x[j], &self.delta_f[j], self.config.relaxation);
                    add_scaled_in_place(&mut x_next, &term, T::ZERO - gamma[j]);
                }
            } else {
                // Degenerate subproblem — clear history and recover with plain relaxation
                self.delta_x.clear();
                self.delta_f.clear();
                if let Some(qr) = &mut self.qr_state {
                    qr.reset();
                }
            }
        }

        self.prev_x = Some(x.clone());
        self.prev_f = Some(f);

        x_next
    }

    /// Solve via normal equations using a pivoted dense solve on the small history matrix.
    fn solve_normal_equations(&self, f: &Array1<T>) -> Option<Array1<T>> {
        let m = self.delta_f.len();
        if m == 0 {
            return None;
        }
        let mut normal = matrix_zeros(m, m);
        let mut rhs = vector_zeros(m);
        for row in 0..m {
            rhs[[row]] = dot(&self.delta_f[row], f);
            for col in 0..m {
                normal[[row, col]] = dot(&self.delta_f[row], &self.delta_f[col]);
            }
        }
        solve_dense_pivoted(normal, rhs, self.config.drop_tolerance)
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

fn solve_dense_pivoted<T: RealField + Copy + FloatElement + std::fmt::Debug>(
    mut matrix: Array2<T>,
    mut rhs: Array1<T>,
    drop_tolerance: T,
) -> Option<Array1<T>> {
    let [rows, cols] = matrix.shape();
    let n = vector_len(&rhs);
    if rows != n || cols != n {
        return None;
    }

    for pivot_col in 0..n {
        let mut pivot_row = pivot_col;
        let mut pivot_abs = NumericElement::abs(matrix[[pivot_col, pivot_col]]);
        for candidate in (pivot_col + 1)..n {
            let candidate_abs = NumericElement::abs(matrix[[candidate, pivot_col]]);
            if candidate_abs > pivot_abs {
                pivot_abs = candidate_abs;
                pivot_row = candidate;
            }
        }
        if pivot_abs < drop_tolerance {
            return None;
        }

        if pivot_row != pivot_col {
            for col in pivot_col..n {
                let pivot_value = matrix[[pivot_col, col]];
                matrix[[pivot_col, col]] = matrix[[pivot_row, col]];
                matrix[[pivot_row, col]] = pivot_value;
            }
            let pivot_rhs = rhs[pivot_col];
            rhs[pivot_col] = rhs[pivot_row];
            rhs[pivot_row] = pivot_rhs;
        }

        for row in (pivot_col + 1)..n {
            let factor = matrix[[row, pivot_col]] / matrix[[pivot_col, pivot_col]];
            matrix[[row, pivot_col]] = T::ZERO;
            for col in (pivot_col + 1)..n {
                matrix[[row, col]] = matrix[[row, col]] - factor * matrix[[pivot_col, col]];
            }
            rhs[[row]] = rhs[[row]] - factor * rhs[[pivot_col]];
        }
    }

    let mut solution = vector_zeros(n);
    for row in (0..n).rev() {
        let diag = matrix[[row, row]];
        if NumericElement::abs(diag) < drop_tolerance {
            return None;
        }
        let mut value = rhs[[row]];
        for col in (row + 1)..n {
            value -= matrix[[row, col]] * solution[[col]];
        }
        solution[[row]] = value / diag;
    }
    Some(solution)
}

#[cfg(test)]
mod tests {
    use super::super::linalg::vector_from_vec;
    use super::*;
    use eunomia::assert_relative_eq;

    fn vec(values: Vec<f64>) -> Array1<f64> {
        vector_from_vec(values)
    }

    fn mat_vec(matrix: &Array2<f64>, vector: &Array1<f64>) -> Array1<f64> {
        let [rows, cols] = matrix.shape();
        assert_eq!(cols, vector_len(vector));
        vector_from_vec(
            (0..rows)
                .map(|row| (0..cols).fold(0.0, |acc, col| acc + matrix[[row, col]] * vector[[col]]))
                .collect(),
        )
    }

    fn make_linear_system() -> (Array2<f64>, Array1<f64>) {
        // G(x) = A x + b, fixed point at x* = [1.0, 2.0]
        let a = Array2::from_shape_vec([2, 2], vec![0.5, 0.1, 0.0, 0.5]).unwrap();
        let b = vec(vec![1.0 - (0.5 * 1.0 + 0.1 * 2.0), 2.0 - 0.5 * 2.0]);
        (a, b)
    }

    fn run_fixed_point(config: AndersonConfig<f64>, iters: usize) -> Array1<f64> {
        let (a, b) = make_linear_system();
        let mut accelerator = AndersonAccelerator::new(config);
        let mut x = vec(vec![0.0, 0.0]);
        for _ in 0..iters {
            let g_x = add_scaled(&mat_vec(&a, &x), &b, 1.0);
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
        assert_relative_eq!(x[[0]], 1.0, epsilon = 1e-8);
        assert_relative_eq!(x[[1]], 2.0, epsilon = 1e-8);
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
        assert_relative_eq!(x_qr[[0]], x_ne[[0]], epsilon = 1e-7);
        assert_relative_eq!(x_qr[[1]], x_ne[[1]], epsilon = 1e-7);
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
        let mut x = vec(vec![0.0, 0.0]);
        for _ in 0..(m + 5) {
            let g_x = add_scaled(&mat_vec(&a, &x), &b, 1.0);
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
        assert!((x[[0]] - 1.0).abs() < 1e-6);
        assert!((x[[1]] - 2.0).abs() < 1e-6);
    }

    /// QR should handle near-collinear history without divergence (ill-conditioning test).
    /// When ΔF columns are nearly identical, the normal equations approach is prone to
    /// numerical rank deficiency — QR should silently discard bad columns.
    #[test]
    fn test_anderson_qr_ill_conditioned_history_graceful() {
        // Very slowly-converging fixed point: contraction ratio 0.99 → lots of similar steps
        let a = Array2::from_shape_vec([2, 2], vec![0.99, 0.0, 0.0, 0.99]).unwrap();
        let b = vec(vec![0.01, 0.02]); // x* = [1.0, 2.0]

        let config = AndersonConfig::<f64> {
            history_depth: 5,
            relaxation: 1.0,
            drop_tolerance: 1e-12,
            method: AndersonMethod::QR,
        };
        let mut acc = AndersonAccelerator::new(config);
        let mut x = vec(vec![0.0, 0.0]);

        // Run until convergence or 200 iterations — must not panic or diverge
        for _ in 0..200 {
            let g_x = add_scaled(&mat_vec(&a, &x), &b, 1.0);
            x = acc.compute_next(&x, &g_x);
            if (x[[0]] - 1.0).abs() < 1e-6 && (x[[1]] - 2.0).abs() < 1e-6 {
                break;
            }
        }
        assert!(
            (x[[0]] - 1.0).abs() < 1e-4,
            "QR Anderson must converge on ill-conditioned problem: x[0]={:.6}",
            x[[0]]
        );
        assert!(
            (x[[1]] - 2.0).abs() < 1e-4,
            "QR Anderson must converge on ill-conditioned problem: x[1]={:.6}",
            x[[1]]
        );
    }

    /// Lockstep invariant (OPEN-033 mitigation): when the QR path silently
    /// rejects a near-collinear column, the parallel `delta_x`/`delta_f`
    /// deques must stay aligned with `qr.q_cols` so subsequent
    /// `compute_next` calls index `(γᵢ, Δxᵢ, Δfᵢ)` triples that refer to
    /// the **same** history entry. A regression here would re-introduce the
    /// desync that the gap_audit flagged.
    #[test]
    fn test_anderson_qr_lockstep_invariant_under_rejection() {
        // Fixed point map that produces ΔF vectors that are highly
        // collinear after iter 0, forcing QR to reject most columns.
        let a = Array2::from_shape_vec([2, 2], vec![0.99, 0.0, 0.0, 0.99]).unwrap();
        // Tiny RHS so the iterates shrink slowly and ΔF stays near-linear.
        let b = vec(vec![1e-12, 2e-12]);

        let config = AndersonConfig::<f64> {
            history_depth: 5,
            relaxation: 1.0,
            // A loose drop tolerance forces *every* new column to be
            // rejected after the first, so we exercise the eviction path
            // many times in a row.
            drop_tolerance: 1e-2,
            method: AndersonMethod::QR,
        };
        let mut acc = AndersonAccelerator::new(config);
        let mut x = vec(vec![1.0, 2.0]);

        for iter in 0..200 {
            let g_x = add_scaled(&mat_vec(&a, &x), &b, 1.0);
            x = acc.compute_next(&x, &g_x);

            // Invariant: every step, the parallel deque lengths must match
            // the QR column count exactly — or be the empty initial state
            // (before the first iteration has produced (Δx, Δf)).
            let qr_cols = acc.qr_state.as_ref().map_or(0, |qr| qr.q_cols.len());
            assert_eq!(
                acc.delta_x.len(),
                qr_cols,
                "iter {}: delta_x.len()={} should equal qr.q_cols.len()={}",
                iter,
                acc.delta_x.len(),
                qr_cols,
            );
            assert_eq!(
                acc.delta_f.len(),
                qr_cols,
                "iter {}: delta_f.len()={} should equal qr.q_cols.len()={}",
                iter,
                acc.delta_f.len(),
                qr_cols,
            );
        }
    }
}
