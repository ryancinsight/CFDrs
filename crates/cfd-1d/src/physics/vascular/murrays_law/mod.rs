//! Murray's Law for optimal vascular bifurcation geometry.
//!
//! ## Theorem: Murray's Law of Minimal Work
//!
//! **Theorem**: The optimal radius r of a vascular segment that minimizes the total
//! biological work W_total required to maintain steady blood flow Q scales as r ∝ Q^(1/3).
//!
//! **Proof Outline**: The total work is the sum of viscous dissipation power (Poiseuille
//! resistance) and the metabolic cost of maintaining the blood volume:
//! ```text
//! W_total = Q² · (8μL / πr⁴) + λ · πr²L
//! ```
//! Setting ∂W/∂r = 0 yields Q ∝ r³.
//!
//! ## Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`law`] | `MurraysLaw<T>` — diameter calculations, deviation, validation |
//! | [`bifurcation`] | `OptimalBifurcation<T>` — symmetric/asymmetric geometry |
//!
//! ## References
//!
//! - Murray, C.D. (1926) *J. Gen. Physiol.* 9, 835-841.
//! - Sherman, T.F. (1981) *J. Gen. Physiol.* 78, 431-453.
//! - Zamir, M. (1978) *J. Theor. Biol.* 74, 227-250.
//! - Revellin, R. et al. (2009) *Theor. Biol. Med. Model.* 6:7.

pub mod bifurcation;
pub mod law;

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use bifurcation::OptimalBifurcation;
pub use law::MurraysLaw;

// ── Non-Newtonian Flow-Split Exponent ───────────────────────────────────────

/// Non-Newtonian flow-split exponent for power-law fluids.
///
/// ## Theorem (Revellin et al. 2009)
///
/// For a power-law fluid with constitutive relation τ = K·γ̇ⁿ, the volumetric
/// flow rate through a circular tube of radius R under pressure gradient ΔP/L is:
///
/// ```text
/// Q = (nπ/(3n+1)) · (ΔP/(2KL))^(1/n) · R^((3n+1)/n)
/// ```
///
/// At a bifurcation with equal pressure drop across both daughter branches,
/// the flow-split ratio is:
///
/// ```text
/// Q₁/Q₂ = (D₁/D₂)^((3n+1)/n)
/// ```
///
/// The exponent m = (3n+1)/n determines how strongly vessel geometry affects
/// the flow distribution in non-Newtonian fluids:
///
/// | n   | Fluid type                | m = (3n+1)/n |
/// |-----|---------------------------|--------------|
/// | 1.0 | Newtonian                 | 4.0          |
/// | 0.9 | Mildly shear-thinning     | 4.11         |
/// | 0.8 | Moderately shear-thinning | 4.25         |
/// | 0.5 | Strongly shear-thinning   | 5.0          |
///
/// **Reference**: Revellin et al. (2009). *Theor. Biol. Med. Model.* 6:7.
///
/// # Panics
/// Panics if `power_law_index_n` is not positive.
pub fn non_newtonian_flow_split_exponent(power_law_index_n: f64) -> f64 {
    assert!(
        power_law_index_n > 0.0,
        "Power-law index must be positive, got {power_law_index_n}"
    );
    (3.0 * power_law_index_n + 1.0) / power_law_index_n
}

#[cfg(test)]
mod tests;
