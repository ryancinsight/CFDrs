//! JFNK fallback for stalled Anderson-accelerated Picard iteration (OPEN-033).
//!
//! # Theorem — Newton-Krylov globalization of non-contractive Picard maps
//!
//! A fixed-point iteration `x_{k+1} = G(x_k)` converges linearly with rate
//! `ρ(G'(x*))` when `ρ < 1` (Banach). When `ρ ≥ 1` along the limit cycle, the
//! Picard iterate oscillates in a bounded envelope without decay and
//! Anderson acceleration — a rank-`m` secant mixing of past iterates — cannot
//! escape the limit cycle because rank-m mixing of cyclic history *is* the
//! cycle. This is the documented failure mode for `NetworkSolver` on
//! ill-conditioned network blueprints with large `R_eff` spread:
//! the *linear* residual `‖Ax − b‖ / ‖b‖` falls to machine ε by iter ~8, but
//! the *fixed-point map* `‖Δx‖ / ‖x‖` hovers at `~1e-3` indefinitely.
//!
//! Newton's iteration `x_{k+1} = x_k − J_F(x_k)^{-1} F(x_k)` with
//! `F(x) := G(x) − x` converges **quadratically** in the basin of any regular
//! root (Kantorovich 1948), and its convergence condition (`det J_F(x*) ≠ 0`)
//! is strictly weaker than the Picard condition (`ρ(G'(x*)) < 1`):
//! a non-contractive Picard map can have a nearby fixed point whose Newton
//! iterate still converges. The Jacobian-free variant (Brown & Saad 1990;
//! Knoll & Keyes 2004) replaces the explicit `J_F` with directional
//! finite-difference products `J v ≈ [F(x + ε v) − F(x)] / ε` solved by an
//! inner Krylov solver, eliminating `O(n²)` storage on large networks.
//!
//! # Fallback trigger: bounded-amplitude stagnation
//!
//! The Picard trajectory enters the JFNK fallback when the recent
//! `solution_change_norm` history shows a non-decaying cycle: amplitude
//! `max/min ≤ stall_ratio` over `stall_window` iterations. This is the
//! Walker & Ni (2011, §3.4) restart criterion restricted to fixed-amplitude
//! cycles — strictly tighter than a residual-only check, since the linear
//! residual stays at machine ε throughout the limit cycle.
//!
//! # Budget discipline
//!
//! `JfnkConfig.max_newton_iterations` and `max_krylov_iterations` are
//! chosen so `warmup + newton × krylov ≤ picard_budget`. The fallback is
//! therefore a *trajectory change* in the exhausted iteration budget, not a
//! budget raise — `tolerance`, `max_iterations`, and the `has_converged_dual`
//! gate are unchanged. `atol := tolerance` makes the JFNK convergence
//! `‖F(x)‖ < tolerance` equivalent in spirit to the Picard dual gate
//! `‖Δx‖ / ‖x‖ < tolerance` once `‖x‖ ≈ O(1)` near the root.
//!
//! # References
//!
//! - Brown, P.N. & Saad, Y. (1990). Hybrid Krylov methods for nonlinear systems
//!   of equations. *SIAM J. Sci. Stat. Comput.* 11(3):450–481.
//! - Eisenstat, S.C. & Walker, H.F. (1996). Choosing the forcing terms in an
//!   inexact Newton method. *SIAM J. Sci. Comput.* 17(1):16–32.
//! - Kantorovich, L.V. (1948). On Newton's method for functional equations.
//!   *Doklady Akademii Nauk SSSR* 59(7):1237–1240.
//! - Knoll, D.A. & Keyes, D.E. (2004). Jacobian-free Newton-Krylov methods:
//!   a survey. *J. Comput. Phys.* 193(2):357–397.
//! - Walker, H.F. & Ni, P. (2011). Anderson acceleration of fixed-point
//!   iteration. *SIAM J. Numer. Anal.* 49(3):809–828.

use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use cfd_math::nonlinear_solver::{JfnkConfig, JfnkSolver};
use eunomia::{FloatElement, NumericElement};
use leto::Array1;

use super::linear_system::{LinearSolverMethod, LinearSystemSolver};
use super::matrix_assembly::MatrixAssembler;
use super::workspace::SolverWorkspace;
use super::{Network, NetworkSolveScalar};
use cfd_core::physics::fluid::FluidTrait;

/// Configuration for the JFNK fallback.
///
/// All defaults preserve the outer `SolverConfig` budget: the fallback is a
/// trajectory change inside the existing iteration budget, never a budget
/// raise or tolerance relaxation.
#[derive(Copy, Clone, Debug)]
pub(super) struct FallbackBudget {
    /// Number of Picard iterations before the stagnation test is consulted.
    /// Lets the limit cycle become well-characterised before judging it
    /// non-decaying.
    pub(super) warmup_iterations: usize,
    /// Window of recent `solution_change_norm` samples used for the
    /// fixed-amplitude test (`max/min ≤ stall_ratio ⟹ stalled`).
    pub(super) stall_window: usize,
    /// Upper bound on the ratio of max-to-min `solution_change_norm` within
    /// `stall_window` that classifies the trajectory as a stalled limit
    /// cycle. Walker & Ni 2011 §3.4 recommend values in `[5, 20]`; 10 is
    /// conservative and well within the band.
    pub(super) stall_ratio: f64,
    /// Maximum outer Newton iterations in the fallback solve. Quadratic
    /// convergence from a Picard-warmed basin needs `~10`; 30 is generous
    /// headroom.
    pub(super) max_newton_iterations: usize,
    /// Maximum inner Krylov (GMRES) iterations per Newton step. Bounded to
    /// keep `newton × krylov` inside the Picard budget.
    pub(super) max_krylov_iterations: usize,
}

impl Default for FallbackBudget {
    fn default() -> Self {
        // With the canonical `SolverConfig { max_iterations = 10000 }`,
        // warmup 200 + newton 30 × krylov 50 = 200 + 1500 = 1700 ≪ 10000.
        // The active test budget roughly halves because each JvP does one
        // extra F-eval, but the production Picard budget is not raised.
        Self {
            warmup_iterations: 200,
            stall_window: 64,
            stall_ratio: 10.0,
            max_newton_iterations: 30,
            max_krylov_iterations: 50,
        }
    }
}

impl FallbackBudget {
    /// Stagnation detector: bounded fixed-amplitude cycle on
    /// `solution_change_norm`.
    ///
    /// Returns `true` once `history` has at least `stall_window` non-trivial
    /// entries whose max-to-min ratio is at most `stall_ratio`. Zero or
    /// non-finite samples are skipped, so a transient NaN or genuinely
    /// converging tail (whose min drops to ε-mach) does not falsely trigger.
    pub(super) fn is_stalling(&self, history: &[f64]) -> bool {
        if history.len() < self.stall_window {
            return false;
        }
        let tail = &history[history.len() - self.stall_window..];
        let mut lo = f64::INFINITY;
        let mut hi = 0.0_f64;
        for &value in tail {
            if value > 0.0 && value.is_finite() {
                lo = lo.min(value);
                hi = hi.max(value);
            }
        }
        if lo == 0.0 || !lo.is_finite() || hi == 0.0 || !hi.is_finite() {
            return false;
        }
        // Walker & Ni 2011 §3.4: a sustained fixed-amplitude envelope with
        // bounded max/min ratio indicates a non-contractive Picard map.
        // The Picard iteration cannot break out of this cycle on its own;
        // a Newton-Krylov globalization reaches the root via `J_F^{-1} F`.
        hi / lo <= self.stall_ratio
    }
}

/// Run JFNK on `F(x) = G(x) − x` starting from the Picard-warmed solution.
///
/// The closure prefers `&Network<T, F>` cloning inside the F-eval (one per
/// JvP) so `network_snapshot` is preserved across calls — JFNK samples
/// `F(x + ε v)` and `F(x)` independently during each Jacobian-vector
/// product. Per-call allocation is acceptable: JFNK only runs after the
/// Picard trajectory has stalled, so this is the cold-path recovery.
///
/// # Errors
///
/// Returns `ConvergenceErrorKind::MaxIterationsExceeded` if JFNK did not
/// drive `‖F‖` below `tolerance` within `max_newton_iterations`.
/// `tolerance` is not modified: pass the same value as the failing Picard
/// gate so the Newton solution satisfies `‖F‖ < tolerance ≈ ‖Δx‖/‖x‖ < tol`.
pub(super) fn jfnk_fallback<T, F>(
    assembler: &MatrixAssembler<T>,
    network_snapshot: &Network<T, F>,
    dirichlet_values: Vec<Option<T>>,
    neumann_sources: Vec<Option<T>>,
    linear_method: LinearSolverMethod,
    tolerance: T,
    fallback: &FallbackBudget,
    warm_solution: &Array1<T>,
) -> Result<Array1<T>>
where
    T: NetworkSolveScalar,
    F: FluidTrait<T> + Clone,
{
    let n = network_snapshot.node_count();
    debug_assert_eq!(
        warm_solution.shape()[0],
        n,
        "JFNK warm-start vector length must match network node count"
    );

    // Pre-allocated workspace and solver reused across all F-evaluations.
    // The boundary-condition classification is constant once the outer
    // loop has assembled it; pass through without re-classifying.
    let mut workspace = SolverWorkspace::new(n);
    workspace.dirichlet_values = dirichlet_values;
    workspace.neumann_sources = neumann_sources;
    let inner_solver = LinearSystemSolver::new()
        .with_method(linear_method)
        .with_tolerance(tolerance)
        .with_max_iterations(fallback.max_krylov_iterations);

    // F(x) := G(x) − x where G(x) is one Picard step.
    //
    // The Picard step is identical to `NetworkSolver`'s own `solve_with_initial_guess`
    // loop: (1) refresh edge flow rates from `x`, (2) assemble the conductance
    // Laplacian, (3) solve `A x = b` for `x_next = G(x)`. Returning
    // `x_next − x` gives the Newton residual of the fixed-point map.
    let picard_residual = |x: &Array1<T>| -> Array1<T> {
        let mut net = network_snapshot.clone();
        net.update_from_solution(x)
            .expect("invariant: network update_from_solution must succeed on JFNK F-eval");
        let matrix = assembler
            .assemble_into(&net, &mut workspace)
            .expect("invariant: assemble_into must succeed on JFNK F-eval");
        let mut x_init = x.clone();
        let g_x = inner_solver
            .solve_with_initial_guess(&matrix, &workspace.rhs, &mut x_init)
            .expect("invariant: linear solve must succeed on JFNK F-eval");
        let mut f = g_x.clone();
        for i in 0..n {
            f[[i]] = g_x[[i]] - x[[i]];
        }
        f
    };

    // atol := tolerance. Newton's `‖F‖ < atol` is structurally equivalent
    // to the Picard `‖Δx‖ / ‖x‖ < tolerance` once the iterate is near `x*`
    // (where `‖x*‖ ≈ O(1)` by network-pressure normalization), so no
    // relaxation of the gate is implied.
    let jfnk_config = JfnkConfig {
        max_newton_iterations: fallback.max_newton_iterations,
        atol: tolerance,
        rtol: tolerance,
        max_krylov_iterations: fallback.max_krylov_iterations,
        krylov_restart: fallback.max_krylov_iterations,
        inner_tol_init: <T as FloatElement>::from_f64(0.5),
        eta_min: <T as FloatElement>::from_f64(1e-4),
        eta_max: <T as FloatElement>::from_f64(0.9),
    };
    let jfnk = JfnkSolver::new(jfnk_config);
    let (solution, convergence) = jfnk.solve(picard_residual, warm_solution.clone())?;
    if !convergence.converged {
        return Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: fallback.max_newton_iterations,
            },
        ));
    }
    Ok(solution)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Empty-history slice is not stalled — fallback should not start.
    #[test]
    fn empty_history_is_not_stalling() {
        let budget = FallbackBudget::default();
        assert!(!budget.is_stalling(&[]));
    }

    /// Brief history (shorter than `stall_window`) is not yet pattern-classified.
    #[test]
    fn short_history_is_not_stalling() {
        let budget = FallbackBudget::default();
        let shorter: Vec<f64> = (0..budget.stall_window - 1).map(|_| 1.0).collect();
        assert!(!budget.is_stalling(&shorter));
    }

    /// Bounded fixed-amplitude cycle is the canonical stalled trajectory.
    #[test]
    fn bounded_amplitude_cycle_is_stalling() {
        let budget = FallbackBudget::default();
        // max/min = 1.0 → ratio 1 ≤ stall_ratio (10) so this is stalled.
        let cycle: Vec<f64> = (0..budget.stall_window)
            .map(|i| if i % 2 == 0 { 1.0e-3 } else { 1.5e-3 })
            .collect();
        assert!(budget.is_stalling(&cycle));
    }

    /// Decaying envelope (typical converging Picard) is NOT stalled.
    #[test]
    fn decaying_envelope_is_not_stalling() {
        let budget = FallbackBudget::default();
        // max/min = 1e0 / 1e-12 = 1e12 ≫ stall_ratio (10) → not stalled.
        let decaying: Vec<f64> = (0..budget.stall_window)
            .map(|i| 10f64.powi(-(i as i32)))
            .collect();
        assert!(!budget.is_stalling(&decaying));
    }

    /// Oscillation at a slightly larger envelope (factor 5x) is still stalled
    /// within Walker-Ni bounds.
    #[test]
    fn narrow_envelope_oscillation_is_stalling() {
        let budget = FallbackBudget::default();
        let narrow: Vec<f64> = (0..budget.stall_window)
            .map(|i| if i % 2 == 0 { 1.0e-3 } else { 5.0e-3 })
            .collect();
        // ratio 5 ≤ 10 → stalled
        assert!(budget.is_stalling(&narrow));
    }

    /// Wide-amplitude oscillation (factor 20x) is treated as productive.
    #[test]
    fn wide_envelope_oscillation_is_not_stalling() {
        let budget = FallbackBudget::default();
        let wide: Vec<f64> = (0..budget.stall_window)
            .map(|i| if i % 2 == 0 { 1.0e-3 } else { 2.0e-2 })
            .collect();
        // ratio 20 > 10 → not stalled
        assert!(!budget.is_stalling(&wide));
    }

    /// Transient non-finite samples are skipped, not interpreted as stagnation.
    #[test]
    fn nan_and_zero_samples_are_skipped() {
        let budget = FallbackBudget::default();
        let mut cycle = vec![0.0_f64; budget.stall_window];
        cycle[0] = f64::NAN;
        cycle[1] = 0.0;
        // All NaN/0 sample → no usable history, not stalled.
        assert!(!budget.is_stalling(&cycle));
    }
}
