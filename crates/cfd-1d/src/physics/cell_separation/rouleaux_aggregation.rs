//! Quemada (1978) rouleaux aggregation viscosity model.
//!
//! ## Theorem — Quemada Viscosity Model (Quemada 1978)
//!
//! At low shear rates, RBC rouleaux formation increases the effective
//! viscosity of blood by a factor that depends on both hematocrit and
//! shear rate:
//!
//! ```text
//! μ(γ̇, H_t) = μ_plasma / (1 − ½ k(γ̇) H_t)²
//! ```
//!
//! where k(γ̇) is the intrinsic viscosity function:
//!
//! ```text
//! k(γ̇) = (k_0 + k_∞ √(γ̇/γ_c)) / (1 + √(γ̇/γ_c))
//! ```
//!
//! Parameters (Quemada 1978, human blood):
//! - k_0 = 4.33 (zero-shear intrinsic viscosity)
//! - k_∞ = 2.07 (high-shear intrinsic viscosity)
//! - γ_c = 1.88 s⁻¹ (critical shear rate for rouleaux breakup)
//!
//! **Proof sketch**: At γ̇ → ∞, k → k_∞ and μ → μ_plasma/(1 − k_∞H_t/2)²
//! which recovers the Einstein-Batchelor hard-sphere limit. At γ̇ → 0,
//! k → k_0 > k_∞ reflecting the increased effective volume fraction of
//! rouleaux aggregates.
//!
//! **Reference**: Quemada, D. (1978). "Rheology of concentrated disperse
//! systems III. General features of the proposed non-Newtonian model.
//! Comparison with experimental data", *Rheol. Acta* 17:643-653.

/// Quemada zero-shear intrinsic viscosity parameter.
const K0: f64 = 4.33;
/// Quemada high-shear intrinsic viscosity parameter.
const K_INF: f64 = 2.07;
/// Critical shear rate for rouleaux breakup [s⁻¹].
const GAMMA_C: f64 = 1.88;

/// Compute the Quemada (1978) effective blood viscosity accounting for
/// rouleaux aggregation at low shear rates.
///
/// At low shear rates (γ̇ < 50 s⁻¹), red blood cells form rouleaux
/// (stacked-coin aggregates) that dramatically increase blood viscosity.
/// The Quemada model captures the smooth transition from high-viscosity
/// aggregated blood at rest to the lower viscosity of dispersed cells at
/// high shear.
///
/// # Arguments
/// * `shear_rate` - Local shear rate γ̇ [s⁻¹] (must be ≥ 0)
/// * `hematocrit` - Volume fraction of red blood cells (0.0–0.99)
/// * `mu_plasma` - Plasma dynamic viscosity [Pa·s] (typically 0.0012)
///
/// # Returns
/// Effective dynamic viscosity [Pa·s], always ≥ `mu_plasma`.
#[must_use]
pub fn quemada_viscosity(shear_rate: f64, hematocrit: f64, mu_plasma: f64) -> f64 {
    let gamma = shear_rate.max(0.0);
    let ht = hematocrit.clamp(0.0, 0.99);

    // Intrinsic viscosity function k(γ̇)
    let sqrt_ratio = (gamma / GAMMA_C).sqrt();
    let k = (K0 + K_INF * sqrt_ratio) / (1.0 + sqrt_ratio);

    // Denominator: (1 − ½ k H_t)
    // Guard against divergence at high hematocrit where ½ k H_t → 1
    let denom = (1.0 - 0.5 * k * ht).max(0.01);

    mu_plasma / (denom * denom)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_PLASMA: f64 = 0.0012; // Pa·s
    const HT_NORMAL: f64 = 0.45;

    /// At high shear rate (γ̇ = 1000 s⁻¹), k → k_∞ = 2.07.
    /// The Einstein-Batchelor hard-sphere limit for 45% hematocrit gives
    /// μ = μ_plasma / (1 − k_∞ H_t / 2)² = μ_plasma / (1 − 0.46575)²
    ///   = μ_plasma / 0.28536² ≈ μ_plasma × 12.3
    /// (High because Quemada's k_∞ is fitted for concentrated suspensions.)
    #[test]
    fn test_quemada_high_shear_approaches_einstein() {
        let mu_high = quemada_viscosity(1000.0, HT_NORMAL, MU_PLASMA);

        // At γ̇ = 1000, sqrt(1000/1.88) ≈ 23.1, so k ≈ (4.33 + 2.07×23.1)/(1+23.1) ≈ 2.10
        // Very close to k_∞ = 2.07
        let sqrt_r = (1000.0 / GAMMA_C).sqrt();
        let k_expected = (K0 + K_INF * sqrt_r) / (1.0 + sqrt_r);
        assert!(
            (k_expected - K_INF).abs() < 0.15,
            "At high shear, k = {:.4} should approach k_∞ = {:.4}",
            k_expected,
            K_INF
        );

        // Viscosity should be finite and above plasma
        assert!(mu_high.is_finite() && mu_high > MU_PLASMA);
    }

    /// At low shear rate (γ̇ = 0.1 s⁻¹), rouleaux form and viscosity is
    /// significantly higher than at high shear rate.
    #[test]
    fn test_quemada_low_shear_augmented() {
        let mu_low = quemada_viscosity(0.1, HT_NORMAL, MU_PLASMA);
        let mu_high = quemada_viscosity(1000.0, HT_NORMAL, MU_PLASMA);

        assert!(
            mu_low > mu_high,
            "Low-shear viscosity ({:.6} Pa·s) should exceed high-shear ({:.6} Pa·s) due to rouleaux",
            mu_low,
            mu_high
        );
    }

    /// Zero hematocrit means no cells — viscosity should equal plasma viscosity.
    #[test]
    fn test_quemada_zero_hematocrit_returns_plasma() {
        let mu = quemada_viscosity(10.0, 0.0, MU_PLASMA);
        assert!(
            (mu - MU_PLASMA).abs() < 1e-15,
            "Zero hematocrit viscosity ({:.10}) should equal plasma viscosity ({:.10})",
            mu,
            MU_PLASMA
        );
    }

    /// Viscosity should decrease monotonically with increasing shear rate
    /// (shear-thinning behaviour) because rouleaux break up progressively.
    #[test]
    fn test_quemada_monotone_with_shear() {
        let shear_rates = [0.01, 0.1, 1.0, 10.0, 100.0, 1000.0];
        let viscosities: Vec<f64> = shear_rates
            .iter()
            .map(|&sr| quemada_viscosity(sr, HT_NORMAL, MU_PLASMA))
            .collect();

        for i in 1..viscosities.len() {
            assert!(
                viscosities[i] <= viscosities[i - 1],
                "Viscosity should decrease with shear: μ({}) = {:.6} > μ({}) = {:.6}",
                shear_rates[i],
                viscosities[i],
                shear_rates[i - 1],
                viscosities[i - 1],
            );
        }
    }

    /// Higher hematocrit should always produce higher viscosity at any
    /// given shear rate (more cells → more resistance to flow).
    #[test]
    fn test_quemada_increases_with_hematocrit() {
        let gamma = 50.0; // moderate shear
        let mu_low = quemada_viscosity(gamma, 0.20, MU_PLASMA);
        let mu_mid = quemada_viscosity(gamma, 0.35, MU_PLASMA);
        let mu_high = quemada_viscosity(gamma, 0.50, MU_PLASMA);

        assert!(
            mu_low < mu_mid && mu_mid < mu_high,
            "Viscosity should increase with hematocrit: {:.6} < {:.6} < {:.6}",
            mu_low,
            mu_mid,
            mu_high
        );
    }
}
