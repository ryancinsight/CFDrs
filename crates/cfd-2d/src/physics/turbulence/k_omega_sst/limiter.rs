//! Menter (2003) SST production limiter.
//!
//! # Theorem — Production Limiter (Menter, Kuntz & Langtry 2003)
//!
//! The turbulence production term in the k-equation is limited to prevent
//! unphysical kinetic energy buildup in stagnation regions:
//!
//! ```text
//! P_k = min(P_k_unlimited, C_lim · β* · ρ · k · ω)
//! ```
//!
//! where `C_lim = 10` (standard value) and `β* = 0.09`.
//!
//! **Proof**: In stagnation regions, the strain rate `S` can be large while
//! the turbulent kinetic energy `k` is small, causing `P_k = μ_t S² >> ε`.
//! The limiter bounds `P_k` to be at most `C_lim` times the dissipation rate,
//! preventing runaway growth that would violate the realizability condition
//! `2k/3 >= max(u'², v'², w'²)`.
//!
//! # References
//! - Menter, F.R., Kuntz, M. & Langtry, R. (2003). "Ten Years of Industrial
//!   Experience with the SST Turbulence Model", *Turbulence, Heat and Mass
//!   Transfer* 4:625-632.

/// Apply the Menter (2003) SST production limiter.
///
/// Returns `min(P_k, C_lim · β* · k · ω)` with `C_lim = 10`.
///
/// # Arguments
/// - `p_k` — unlimited production term `μ_t · 2 · S_ij · S_ij`
/// - `k` — turbulent kinetic energy [m²/s²]
/// - `omega` — specific dissipation rate [1/s]
/// - `beta_star` — model constant (0.09 for SST)
#[inline]
#[must_use]
pub fn limit_production(p_k: f64, k: f64, omega: f64, beta_star: f64) -> f64 {
    let c_lim = 10.0;
    p_k.min(c_lim * beta_star * k * omega)
}
