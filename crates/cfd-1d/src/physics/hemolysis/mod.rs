//! Haemolysis index models for millifluidic and microfluidic flow.
//!
//! Provides the Giersiepen (1990) power-law haemolysis index model and a
//! conservative cavitation-amplification correction for SDT (Sonodynamic
//! Therapy) millifluidic devices.
//!
//! # Giersiepen (1990) Model
//!
//! The empirical power-law model relates haemoglobin release to shear stress
//! and exposure duration:
//!
//! ```text
//!   HI = C · t^α · τ^β
//! ```
//!
//! where:
//! - `C  = 3.62 × 10⁻⁵`  (fit constant)
//! - `α  = 0.765`         (time exponent)
//! - `β  = 1.991`         (shear exponent)
//! - `t` = exposure duration [s]
//! - `τ` = wall shear stress [Pa]
//!
//! Reference: Giersiepen M. et al. (1990) *Estimation of shear stress-related
//! blood damage in heart valve prostheses*. Int. J. Artif. Organs 13(5):300–306.
//!
//! # Cavitation Amplification
//!
//! In SDT millifluidic devices, acoustic bubble collapse generates micro-jets
//! and shockwaves that cause RBC membrane damage independently of the
//! macroscopic shear stress.  A conservative 3× amplification factor at full
//! cavitation potential is applied:
//!
//! ```text
//!   HI_amplified = HI_base × (1 + 3 × cav_potential)
//! ```

pub mod acoustic_radiation;

// ── Giersiepen model constants (re-exported from cfd-core SSOT) ───────────────
//
// These are the single source of truth from cfd-core, aliased here for backward
// compatibility with existing callers.  All three fidelity levels (1D/2D/3D)
// now share identical constants, ensuring cross-fidelity HI comparisons are valid.

/// Giersiepen (1990) fit constant C — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_C`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_C as GIERSIEPEN_C;

/// Giersiepen (1990) time exponent α — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_TIME`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_TIME as GIERSIEPEN_ALPHA;

/// Giersiepen (1990) shear exponent β — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_STRESS`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_STRESS as GIERSIEPEN_BETA;

/// Conservative cavitation amplification slope — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::CAVITATION_HI_SLOPE`].
pub use cfd_core::physics::hemolysis::CAVITATION_HI_SLOPE;

// ── Public API ────────────────────────────────────────────────────────────────

/// Giersiepen (1990) haemolysis index: `HI = C · t^α · τ^β`.
///
/// Thin wrapper delegating to [`cfd_core::physics::hemolysis::HemolysisModel::giersiepen_millifluidic`]
/// so that the single source of truth for all model constants lives in cfd-core.
/// Returns `0.0` for non-positive inputs (no damage when shear or time is
/// zero or negative) and `0.0` on conversion error (never expected in practice).
///
/// # Arguments
///
/// * `shear_pa`   — wall shear stress [Pa]
/// * `duration_s` — exposure duration [s]
///
/// # Example
///
/// ```
/// use cfd_1d::giersiepen_hi;
/// let hi = giersiepen_hi(50.0, 0.01);
/// assert!((hi - 2.578_444_864_181_722_3e-3).abs() < 1e-15);
/// ```
#[inline]
#[must_use]
pub fn giersiepen_hi(shear_pa: f64, duration_s: f64) -> f64 {
    cfd_core::physics::hemolysis::HemolysisModel::giersiepen_millifluidic()
        .damage_index(shear_pa, duration_s)
        .unwrap_or(0.0)
}

/// Amplify a baseline Giersiepen HI by the local cavitation potential.
///
/// Thin wrapper delegating to [`cfd_core::physics::hemolysis::HemolysisModel::cavitation_amplified`]
/// so that the single source of truth for the amplification slope lives in cfd-core.
///
/// Bubble collapse generates micro-jets (localised shear amplification) and
/// pressure shockwaves that damage RBC membranes independently of macroscopic
/// steady shear.
///
/// Conservative model: 3× amplification at full cavitation potential
/// (`cav_potential = 1.0`).
///
/// ```text
///   HI_amplified = base_hi × (1 + 3 × cav_potential.clamp(0, 1))
/// ```
///
/// # Arguments
///
/// * `base_hi`          — baseline Giersiepen HI (from [`giersiepen_hi`])
/// * `cav_potential`    — cavitation potential ∈ [0, 1]; 0 = no cavitation
///
/// # Example
///
/// ```
/// use cfd_1d::{cavitation_amplified_hi, giersiepen_hi};
/// let base = giersiepen_hi(100.0, 1.0);
/// let amplified = cavitation_amplified_hi(base, 1.0);
/// assert!((amplified - base * 4.0).abs() < 1e-15);  // 1 + 3×1.0 = 4×
/// ```
#[inline]
#[must_use]
pub fn cavitation_amplified_hi(base_hi: f64, cav_potential: f64) -> f64 {
    cfd_core::physics::hemolysis::HemolysisModel::cavitation_amplified(base_hi, cav_potential)
}

// ── Taskin (2012) strain-based hemolysis model ──────────────────────────────

/// Taskin (2012) calibration constant `C_T = 1.228 × 10⁻⁵`.
pub const TASKIN_C: f64 = 1.228e-5;

/// Taskin (2012) shear stress exponent `β_T = 1.9918`.
///
/// Slightly different from the Giersiepen value of 1.991, reflecting the
/// re-calibration against Lagrangian particle tracking data.
pub const TASKIN_BETA: f64 = 1.9918;

/// Taskin (2012) strain-based haemolysis index (single-segment evaluation).
///
/// Computes the haemolysis index using the Taskin integral-form model for a
/// single channel segment with constant shear stress:
///
/// ```text
/// HI_Taskin = C_T · τ^β_T · t
/// ```
///
/// This is the single-segment (constant shear) evaluation of the full
/// integral form:
///
/// ```text
/// HI_Taskin = C_T · ∫ τ(t)^β_T dt
/// ```
///
/// ## Theorem: Strain-Rate Path Dependence (Taskin et al. 2012)
///
/// The Taskin model captures cumulative damage along the flow path by integrating
/// the instantaneous shear-stress contribution over time, rather than using a
/// single power-law evaluation as in Giersiepen (1990). This integral form
/// correctly accounts for varying shear exposure histories and yields more
/// accurate hemolysis predictions for complex flow geometries where shear stress
/// varies along particle trajectories.
///
/// For constant shear stress, the integral reduces to `C_T · τ^β_T · t`, which
/// differs from the Giersiepen model (`C · t^α · τ^β`) in that the time
/// dependence is linear rather than sub-linear (`α = 0.765` in Giersiepen).
///
/// ## Reference
///
/// Taskin, M. E. et al. (2012). Evaluation of Eulerian and Lagrangian Models
/// for Hemolysis Estimation. *ASAIO J.*, 58(4), 363–372.
///
/// # Arguments
///
/// * `shear_stress` — wall shear stress [Pa]
/// * `exposure_time` — exposure duration [s]
///
/// # Returns
///
/// Haemolysis index (dimensionless). Returns `0.0` for non-positive inputs.
#[inline]
#[must_use]
pub fn taskin_hi(shear_stress: f64, exposure_time: f64) -> f64 {
    if shear_stress <= 0.0 || exposure_time <= 0.0 {
        return 0.0;
    }
    TASKIN_C * shear_stress.abs().max(1e-30).powf(TASKIN_BETA) * exposure_time
}

// ── Convenience struct ────────────────────────────────────────────────────────

/// Composite haemolysis exposure combining shear stress, duration, and
/// cavitation potential.
///
/// Use [`HemolysisExposure::compute_index`] to obtain the amplified HI for
/// a single exposure event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HemolysisExposure {
    /// Wall shear stress at the exposure site [Pa].
    pub shear_pa: f64,
    /// Duration of the exposure event [s].
    pub duration_s: f64,
    /// Local cavitation potential ∈ [0, 1]; 0 for non-venturi regions.
    pub cavitation_potential: f64,
}

impl HemolysisExposure {
    /// Create a new exposure event.
    #[must_use]
    pub fn new(shear_pa: f64, duration_s: f64, cavitation_potential: f64) -> Self {
        Self {
            shear_pa,
            duration_s,
            cavitation_potential,
        }
    }

    /// Create an exposure event with no cavitation contribution.
    #[must_use]
    pub fn shear_only(shear_pa: f64, duration_s: f64) -> Self {
        Self::new(shear_pa, duration_s, 0.0)
    }

    /// Compute the amplified haemolysis index for this exposure event.
    ///
    /// Combines [`giersiepen_hi`] with [`cavitation_amplified_hi`].
    #[must_use]
    pub fn compute_index(&self) -> f64 {
        cavitation_amplified_hi(
            giersiepen_hi(self.shear_pa, self.duration_s),
            self.cavitation_potential,
        )
    }
}

// ── Sonosensitizer Activation Kinetics (Rosenthal 2004) ──────────────────────

/// Activation rate constant for Hematoporphyrin \[s⁻¹\] (Rosenthal 2004).
pub const SENSITIZER_K_ACT_HEMATOPORPHYRIN: f64 = 0.5;

/// Activation rate constant for Chlorin e6 \[s⁻¹\].
pub const SENSITIZER_K_ACT_CHLORIN_E6: f64 = 0.8;

/// Sonosensitizer activation efficiency as a function of transit time.
///
/// ## Theorem — First-Order Activation Kinetics (Rosenthal 2004)
///
/// The fraction of sonosensitizer activated during transit through a
/// cavitation zone is governed by first-order kinetics:
///
/// ```text
/// η_act = 1 − exp(−k_act · I_cav · t_transit)
/// ```
///
/// where:
/// - k_act ≈ 0.1–1.0 s⁻¹ (activation rate constant, drug-dependent)
/// - I_cav = cavitation intensity (dimensionless, 0–1)
/// - t_transit = residence time in cavitation zone \[s\]
///
/// For short transits (k·I·t << 1): η_act ≈ k·I·t (linear regime)
/// For long transits (k·I·t >> 1): η_act → 1 (saturation)
///
/// **Reference**: Rosenthal, I., Sostaric, J.Z. & Riesz, P. (2004).
/// "Sonodynamic therapy – a review of the synergistic effects of drugs
/// and ultrasound", *Ultrason. Sonochem.* 11:349-363.
///
/// # Arguments
///
/// * `k_act` — activation rate constant \[s⁻¹\]
/// * `cavitation_intensity` — I_cav (0–1)
/// * `transit_time_s` — residence time in cavitation zone \[s\]
///
/// # Returns
///
/// Activation fraction η_act ∈ \[0, 1\].
///
/// # Example
///
/// ```
/// use cfd_1d::sonosensitizer_activation_efficiency;
/// let eta = sonosensitizer_activation_efficiency(0.5, 0.8, 2.0);
/// assert!(eta > 0.0 && eta < 1.0);
/// ```
#[inline]
#[must_use]
pub fn sonosensitizer_activation_efficiency(
    k_act: f64,
    cavitation_intensity: f64,
    transit_time_s: f64,
) -> f64 {
    (1.0 - (-k_act * cavitation_intensity * transit_time_s).max(-500.0).exp()).clamp(0.0, 1.0)
}

// ── Rayleigh-Plesset Cavitation Bubble Models (Rayleigh 1917) ────────────────

/// Empirical collapse amplification coefficient α ≈ 3.0.
pub const RAYLEIGH_ALPHA: f64 = 3.0;

/// Standard atmospheric pressure \[Pa\] used as reference for amplification.
pub const P_REF_ATMOSPHERIC: f64 = 101_325.0;

/// Rayleigh inertial collapse time for a cavitation bubble \[s\].
///
/// ## Theorem — Rayleigh Collapse (Rayleigh 1917)
///
/// For a spherical bubble collapsing from radius R_max under external
/// pressure p_inf (neglecting gas content and surface tension):
///
/// ```text
/// t_collapse = 0.915 · R_max · √(ρ / p_inf)
/// ```
///
/// **Reference**: Rayleigh, Lord (1917). "On the pressure developed in a
/// liquid during the collapse of a spherical cavity",
/// *Phil. Mag.* 34:94-98.
///
/// # Arguments
///
/// * `r_max` — maximum bubble radius \[m\]
/// * `rho` — liquid density \[kg/m³\]
/// * `p_inf` — external (driving) pressure \[Pa\]
///
/// # Example
///
/// ```
/// use cfd_1d::rayleigh_collapse_time;
/// let t = rayleigh_collapse_time(50e-6, 1060.0, 101325.0);
/// assert!(t > 0.0);
/// ```
#[inline]
#[must_use]
pub fn rayleigh_collapse_time(r_max: f64, rho: f64, p_inf: f64) -> f64 {
    if r_max <= 0.0 || rho <= 0.0 || p_inf <= 0.0 {
        return 0.0;
    }

    0.915 * r_max * (rho / p_inf).sqrt()
}

/// Micro-jet velocity from Rayleigh bubble collapse \[m/s\].
///
/// The collapse generates a micro-jet with velocity:
///
/// ```text
/// v_jet ≈ √(2 · p_inf / ρ)
/// ```
///
/// # Arguments
///
/// * `p_inf` — external pressure \[Pa\]
/// * `rho` — liquid density \[kg/m³\]
///
/// # Example
///
/// ```
/// use cfd_1d::collapse_jet_velocity;
/// let v = collapse_jet_velocity(101325.0, 1060.0);
/// assert!(v > 10.0, "Jet velocity at atmospheric pressure should exceed 10 m/s");
/// ```
#[inline]
#[must_use]
pub fn collapse_jet_velocity(p_inf: f64, rho: f64) -> f64 {
    if p_inf <= 0.0 || rho <= 0.0 {
        return 0.0;
    }

    (2.0 * p_inf / rho).sqrt()
}

/// Hemolysis amplification factor from cavitation bubble collapse.
///
/// The amplification scales with collapse intensity:
///
/// ```text
/// A_collapse = 1 + α · (R_max / R_0)² · (p_inf / p_ref)
/// ```
///
/// where α ≈ 3.0 (empirical), R_0 is equilibrium radius,
/// p_ref = 101 325 Pa (1 atm).
///
/// # Arguments
///
/// * `r_max` — maximum bubble radius \[m\]
/// * `r_0` — equilibrium bubble radius \[m\]
/// * `p_inf` — external pressure \[Pa\]
///
/// # Example
///
/// ```
/// use cfd_1d::cavitation_hemolysis_amplification;
/// let a = cavitation_hemolysis_amplification(100e-6, 50e-6, 101325.0);
/// assert!(a > 1.0, "Amplification should exceed unity");
/// ```
#[inline]
#[must_use]
pub fn cavitation_hemolysis_amplification(r_max: f64, r_0: f64, p_inf: f64) -> f64 {
    if r_max <= 0.0 || r_0 <= 0.0 || p_inf <= 0.0 {
        return 1.0;
    }

    1.0 + RAYLEIGH_ALPHA * (r_max / r_0).powi(2) * (p_inf / P_REF_ATMOSPHERIC)
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── giersiepen_hi ────────────────────────────────────────────────────────

    #[test]
    fn giersiepen_zero_shear_returns_zero() {
        assert_eq!(giersiepen_hi(0.0, 1.0), 0.0);
    }

    #[test]
    fn giersiepen_zero_duration_returns_zero() {
        assert_eq!(giersiepen_hi(100.0, 0.0), 0.0);
    }

    #[test]
    fn giersiepen_negative_inputs_return_zero() {
        assert_eq!(giersiepen_hi(-10.0, 1.0), 0.0);
        assert_eq!(giersiepen_hi(10.0, -1.0), 0.0);
    }

    #[test]
    fn giersiepen_reference_value_at_100pa_1s() {
        // HI = 3.62e-5 × 1^0.765 × 100^1.991
        let expected =
            GIERSIEPEN_C * 1.0_f64.powf(GIERSIEPEN_ALPHA) * 100.0_f64.powf(GIERSIEPEN_BETA);
        let hi = giersiepen_hi(100.0, 1.0);
        assert!(
            (hi - expected).abs() < 1e-15,
            "got {hi}, expected {expected}"
        );
    }

    #[test]
    fn giersiepen_increases_with_shear() {
        assert!(giersiepen_hi(200.0, 1.0) > giersiepen_hi(100.0, 1.0));
    }

    #[test]
    fn giersiepen_increases_with_duration() {
        assert!(giersiepen_hi(100.0, 2.0) > giersiepen_hi(100.0, 1.0));
    }

    // ── cavitation_amplified_hi ──────────────────────────────────────────────

    #[test]
    fn cavitation_zero_potential_no_amplification() {
        let base = 0.5;
        assert!((cavitation_amplified_hi(base, 0.0) - base).abs() < 1e-15);
    }

    #[test]
    fn cavitation_full_potential_four_times_base() {
        // (1 + 3 × 1.0) = 4.0
        let base = 0.1;
        assert!((cavitation_amplified_hi(base, 1.0) - base * 4.0).abs() < 1e-15);
    }

    #[test]
    fn cavitation_clamped_above_one() {
        let base = 0.1;
        // cav_potential > 1 should be clamped to 1
        assert!(
            (cavitation_amplified_hi(base, 2.0) - cavitation_amplified_hi(base, 1.0)).abs() < 1e-15
        );
    }

    #[test]
    fn cavitation_clamped_below_zero() {
        let base = 0.1;
        // cav_potential < 0 should be clamped to 0
        assert!((cavitation_amplified_hi(base, -1.0) - base).abs() < 1e-15);
    }

    // ── HemolysisExposure ────────────────────────────────────────────────────

    #[test]
    fn exposure_shear_only_matches_giersiepen() {
        let exposure = HemolysisExposure::shear_only(80.0, 0.5);
        let expected = giersiepen_hi(80.0, 0.5);
        assert!((exposure.compute_index() - expected).abs() < 1e-15);
    }

    #[test]
    fn exposure_with_cavitation_exceeds_base() {
        let base = HemolysisExposure::shear_only(80.0, 0.5).compute_index();
        let cav = HemolysisExposure::new(80.0, 0.5, 0.5).compute_index();
        assert!(cav > base);
    }

    // ── taskin_hi ────────────────────────────────────────────────────────────

    #[test]
    fn test_taskin_zero_inputs() {
        assert_eq!(taskin_hi(0.0, 1.0), 0.0);
        assert_eq!(taskin_hi(100.0, 0.0), 0.0);
        assert_eq!(taskin_hi(-10.0, 1.0), 0.0);
        assert_eq!(taskin_hi(10.0, -1.0), 0.0);
    }

    #[test]
    fn test_taskin_reference_value() {
        // HI = C_T × τ^β_T × t = 1.228e-5 × 100^1.9918 × 1.0
        let expected = TASKIN_C * 100.0_f64.powf(TASKIN_BETA) * 1.0;
        let hi = taskin_hi(100.0, 1.0);
        assert!(
            (hi - expected).abs() < 1e-15,
            "got {hi}, expected {expected}"
        );
        // Sanity: should be on the order of 0.1 (100 Pa is significant shear)
        assert!(hi > 0.01 && hi < 10.0, "unexpected magnitude: {hi}");
    }

    #[test]
    fn test_taskin_monotonic() {
        // HI increases with shear stress
        assert!(taskin_hi(200.0, 1.0) > taskin_hi(100.0, 1.0));
        // HI increases with exposure time
        assert!(taskin_hi(100.0, 2.0) > taskin_hi(100.0, 1.0));
    }

    #[test]
    fn test_taskin_vs_giersiepen() {
        // At the same inputs, Taskin and Giersiepen should give different values
        // because they use different constants and time dependence.
        let tau = 100.0;
        let t = 1.0;
        let hi_g = giersiepen_hi(tau, t);
        let hi_t = taskin_hi(tau, t);
        // Both should be positive
        assert!(hi_g > 0.0);
        assert!(hi_t > 0.0);
        // They should differ (different model calibration)
        assert!(
            (hi_g - hi_t).abs() > 1e-10,
            "Giersiepen ({hi_g}) and Taskin ({hi_t}) should give different values"
        );
    }

    // ── sonosensitizer_activation_efficiency ────────────────────────────────

    #[test]
    fn test_activation_zero_time() {
        let eta = sonosensitizer_activation_efficiency(0.5, 0.8, 0.0);
        assert!(
            eta.abs() < 1e-15,
            "Zero transit time should give zero activation, got {eta}"
        );
    }

    #[test]
    fn test_activation_saturates_long_transit() {
        let eta = sonosensitizer_activation_efficiency(0.5, 1.0, 100.0);
        assert!(
            (eta - 1.0).abs() < 1e-10,
            "Very long transit should saturate to 1.0, got {eta}"
        );
    }

    #[test]
    fn test_activation_linear_short_transit() {
        // For k·I·t = 0.01 (small), η ≈ k·I·t (first-order Taylor)
        let k = 0.5;
        let i_cav = 0.1;
        let t = 0.2; // k·I·t = 0.01
        let eta = sonosensitizer_activation_efficiency(k, i_cav, t);
        let linear_approx = k * i_cav * t;
        assert!(
            (eta - linear_approx).abs() < 1e-4,
            "For small k·I·t, η ({eta}) should approximate k·I·t ({linear_approx})"
        );
    }

    #[test]
    fn test_activation_increases_with_intensity() {
        let eta_low = sonosensitizer_activation_efficiency(0.5, 0.2, 1.0);
        let eta_high = sonosensitizer_activation_efficiency(0.5, 0.8, 1.0);
        assert!(
            eta_high > eta_low,
            "Higher intensity should give higher activation: {eta_high} vs {eta_low}"
        );
    }

    #[test]
    fn test_activation_known_constants() {
        // Hematoporphyrin rate constant
        assert!(
            (SENSITIZER_K_ACT_HEMATOPORPHYRIN - 0.5).abs() < 1e-15,
            "Hematoporphyrin k_act should be 0.5"
        );
        // Chlorin e6 rate constant
        assert!(
            (SENSITIZER_K_ACT_CHLORIN_E6 - 0.8).abs() < 1e-15,
            "Chlorin e6 k_act should be 0.8"
        );
    }

    // ── rayleigh_collapse_time ──────────────────────────────────────────────

    #[test]
    fn test_rayleigh_collapse_time_positive() {
        let t = rayleigh_collapse_time(50e-6, 1060.0, 101_325.0);
        assert!(t > 0.0, "Collapse time must be positive, got {t}");
    }

    #[test]
    fn test_rayleigh_collapse_time_scales_with_radius() {
        let t1 = rayleigh_collapse_time(50e-6, 1060.0, 101_325.0);
        let t2 = rayleigh_collapse_time(100e-6, 1060.0, 101_325.0);
        // t ∝ R_max → doubling radius doubles collapse time
        let ratio = t2 / t1;
        assert!(
            (ratio - 2.0).abs() < 1e-10,
            "Doubling R_max should double collapse time, got ratio {ratio}"
        );
    }

    #[test]
    fn test_rayleigh_collapse_time_order_of_magnitude() {
        // For R_max = 50 µm in blood at atmospheric pressure:
        // t = 0.915 × 50e-6 × √(1060/101325) ≈ 4.68e-6 s (~5 µs)
        let t = rayleigh_collapse_time(50e-6, 1060.0, 101_325.0);
        assert!(
            t > 1e-7 && t < 1e-4,
            "Collapse time for 50 µm bubble should be on microsecond scale, got {t}"
        );
    }

    #[test]
    fn test_rayleigh_and_jet_invalid_inputs_return_zero() {
        assert_eq!(rayleigh_collapse_time(0.0, 1060.0, 101_325.0), 0.0);
        assert_eq!(rayleigh_collapse_time(50e-6, 0.0, 101_325.0), 0.0);
        assert_eq!(rayleigh_collapse_time(50e-6, 1060.0, 0.0), 0.0);
        assert_eq!(collapse_jet_velocity(0.0, 1060.0), 0.0);
        assert_eq!(collapse_jet_velocity(101_325.0, 0.0), 0.0);
    }

    // ── collapse_jet_velocity ───────────────────────────────────────────────

    #[test]
    fn test_jet_velocity_at_atmospheric() {
        let v = collapse_jet_velocity(101_325.0, 1060.0);
        // v = √(2 × 101325 / 1060) ≈ 13.83 m/s
        let expected = (2.0 * 101_325.0 / 1060.0_f64).sqrt();
        assert!(
            (v - expected).abs() < 1e-10,
            "Jet velocity should be ~13.8 m/s, got {v}"
        );
        assert!(
            (v - 13.83).abs() < 0.1,
            "Jet velocity at atmospheric pressure should be ~13.8 m/s, got {v}"
        );
    }

    #[test]
    fn test_jet_velocity_increases_with_pressure() {
        let v1 = collapse_jet_velocity(101_325.0, 1060.0);
        let v2 = collapse_jet_velocity(2.0 * 101_325.0, 1060.0);
        assert!(
            v2 > v1,
            "Higher pressure should give higher jet velocity"
        );
    }

    // ── cavitation_hemolysis_amplification ──────────────────────────────────

    #[test]
    fn test_amplification_unity_at_equilibrium() {
        // When R_max = R_0 and p_inf = p_ref → A = 1 + α × 1 × 1 = 1 + 3 = 4
        let a = cavitation_hemolysis_amplification(50e-6, 50e-6, P_REF_ATMOSPHERIC);
        let expected = 1.0 + RAYLEIGH_ALPHA;
        assert!(
            (a - expected).abs() < 1e-10,
            "At equilibrium with atmospheric pressure, A should be {expected}, got {a}"
        );
    }

    #[test]
    fn test_amplification_increases_with_bubble_size() {
        let a1 = cavitation_hemolysis_amplification(50e-6, 50e-6, 101_325.0);
        let a2 = cavitation_hemolysis_amplification(100e-6, 50e-6, 101_325.0);
        assert!(
            a2 > a1,
            "Larger bubble (R_max) should produce stronger amplification: {a2} vs {a1}"
        );
    }

    #[test]
    fn test_amplification_always_above_unity() {
        // Even for small bubbles, amplification should exceed 1.0
        let a = cavitation_hemolysis_amplification(1e-6, 1e-6, 50_000.0);
        assert!(
            a > 1.0,
            "Amplification factor must always exceed 1.0, got {a}"
        );
    }

    #[test]
    fn test_amplification_invalid_inputs_return_unity() {
        assert_eq!(cavitation_hemolysis_amplification(0.0, 1e-6, 101_325.0), 1.0);
        assert_eq!(cavitation_hemolysis_amplification(1e-6, 0.0, 101_325.0), 1.0);
        assert_eq!(cavitation_hemolysis_amplification(1e-6, 1e-6, 0.0), 1.0);
    }

    #[test]
    fn test_amplification_scales_with_r_squared() {
        let r0 = 25e-6;
        let a1 = cavitation_hemolysis_amplification(50e-6, r0, 101_325.0);
        let a2 = cavitation_hemolysis_amplification(100e-6, r0, 101_325.0);
        // A = 1 + α·(R/R0)²·(p/p_ref)
        // Doubling R_max → (R/R0)² quadruples → the α·(R/R0)²·(p/p_ref) term quadruples
        let term1 = a1 - 1.0;
        let term2 = a2 - 1.0;
        let ratio = term2 / term1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "Doubling R_max should quadruple the amplification term, got ratio {ratio}"
        );
    }
}
