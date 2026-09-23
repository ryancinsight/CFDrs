// ── Constraint helpers ────────────────────────────────────────────────────────

/// Smooth sigmoid penalty based on a signed feasibility margin.
///
/// Returns `1.0` when the margin is `≥ +0.1` (well inside the feasible region),
/// `0.5` at the boundary (`margin = 0`), and `0.0` when margin `≤ −0.1`
/// (deeply infeasible).
///
/// # Theorem
/// `sigmoid_penalty(m)` is bounded and monotone:
/// 1. `0 ≤ sigmoid_penalty(m) ≤ 1` for all real `m`
/// 2. if `m1 ≤ m2`, then `sigmoid_penalty(m1) ≤ sigmoid_penalty(m2)`
///
/// **Proof sketch**
/// The function is affine (`0.5 + 5m`) followed by `clamp(0, 1)`.
/// Clamping maps all inputs into `\[0,1]` and preserves monotonicity because
/// both the affine map and clamp are monotone non-decreasing.
#[inline]
#[must_use]
pub fn sigmoid_penalty(margin: f64) -> f64 {
    (0.5 + margin * 5.0).clamp(0.0, 1.0)
}
