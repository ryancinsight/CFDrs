//! Spalding's (1961) implicit universal law of the wall.
//!
//! ## Theorem — Universal Law of the Wall (Spalding 1961)
//!
//! A single implicit relation valid across ALL y⁺ ranges (viscous sublayer,
//! buffer layer, and log-law region):
//!
//! ```text
//! y⁺ = u⁺ + exp(−κB) · [exp(κu⁺) − 1 − κu⁺ − (κu⁺)²/2 − (κu⁺)³/6]
//! ```
//!
//! where κ = 0.41 (von Kármán constant) and B = 5.0 (log-law intercept).
//!
//! This formula provides C¹ continuity (no derivative discontinuities at
//! y⁺ = 5 or y⁺ = 30) and is self-consistent for implicit wall BC solvers.
//!
//! **Proof**: The formula reduces to u⁺ = y⁺ for y⁺ → 0 (Taylor expansion
//! of exp(κu⁺) cancels the leading-order terms) and to u⁺ = ln(y⁺)/κ + B
//! for y⁺ → ∞ (exponential term dominates).
//!
//! ## Reference
//!
//! Spalding, D.B. (1961). "A Single Formula for the Law of the Wall",
//! *J. Appl. Mech.* 28(3):455–458.

/// Given y⁺, returns u⁺ by Newton-Raphson iteration of Spalding's formula.
pub fn spalding_u_plus(y_plus: f64) -> f64 {
    const KAPPA_S: f64 = 0.41;
    const B_S: f64 = 5.0;
    const EXP_NEG_KB: f64 = 0.1108031584; // exp(-0.41 * 5.0) precomputed
    const TOL: f64 = 1e-10;
    const MAX_ITER: usize = 50;

    // Initial guess
    let mut u_plus = if y_plus < 5.0 {
        y_plus
    } else {
        y_plus.ln() / KAPPA_S + B_S
    };

    for _ in 0..MAX_ITER {
        let ku = KAPPA_S * u_plus;
        let exp_ku = ku.exp();

        // F(u⁺) = u⁺ + exp(-κB)·[exp(κu⁺) - 1 - κu⁺ - (κu⁺)²/2 - (κu⁺)³/6] - y⁺
        let f =
            u_plus + EXP_NEG_KB * (exp_ku - 1.0 - ku - ku * ku / 2.0 - ku * ku * ku / 6.0) - y_plus;

        // F'(u⁺) = 1 + exp(-κB)·κ·[exp(κu⁺) - 1 - κu⁺ - (κu⁺)²/2]
        let df = 1.0 + EXP_NEG_KB * KAPPA_S * (exp_ku - 1.0 - ku - ku * ku / 2.0);

        if df.abs() < 1e-30 {
            break;
        }

        let delta = f / df;
        u_plus -= delta;

        if f.abs() < TOL {
            break;
        }
    }

    u_plus
}
