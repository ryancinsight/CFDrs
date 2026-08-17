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
//! - `t` = exposure duration \[s]
//! - `τ` = wall shear stress \[Pa]
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
mod dynamics;
mod models;

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
pub use dynamics::{
    cavitation_hemolysis_amplification, collapse_jet_velocity, rayleigh_collapse_time,
    HemolysisExposure, P_REF_ATMOSPHERIC, RAYLEIGH_ALPHA,
};
pub use models::{
    cavitation_amplified_hi, giersiepen_hi, sonosensitizer_activation_efficiency, taskin_hi,
    SENSITIZER_K_ACT_CHLORIN_E6, SENSITIZER_K_ACT_HEMATOPORPHYRIN, TASKIN_BETA, TASKIN_C,
};

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::quantities::{Pressure, Time};

    fn shear(value: f64) -> Pressure {
        Pressure::from_base(value)
    }

    fn duration(value: f64) -> Time {
        Time::from_base(value)
    }

    // ── giersiepen_hi ────────────────────────────────────────────────────────

    #[test]
    fn giersiepen_zero_shear_returns_zero() {
        assert_eq!(giersiepen_hi(shear(0.0), duration(1.0)), 0.0);
    }

    #[test]
    fn giersiepen_zero_duration_returns_zero() {
        assert_eq!(giersiepen_hi(shear(100.0), duration(0.0)), 0.0);
    }

    #[test]
    fn giersiepen_negative_inputs_return_zero() {
        assert_eq!(giersiepen_hi(shear(-10.0), duration(1.0)), 0.0);
        assert_eq!(giersiepen_hi(shear(10.0), duration(-1.0)), 0.0);
    }

    #[test]
    fn giersiepen_nan_inputs_propagate_nan() {
        assert!(giersiepen_hi(shear(f64::NAN), duration(1.0)).is_nan());
        assert!(giersiepen_hi(shear(1.0), duration(f64::NAN)).is_nan());
    }

    #[test]
    fn giersiepen_reference_value_at_100pa_1s() {
        // HI = 3.62e-5 × 1^0.765 × 100^1.991
        let expected =
            GIERSIEPEN_C * 1.0_f64.powf(GIERSIEPEN_ALPHA) * 100.0_f64.powf(GIERSIEPEN_BETA);
        let hi = giersiepen_hi(shear(100.0), duration(1.0));
        assert!(
            (hi - expected).abs() < 1e-15,
            "got {hi}, expected {expected}"
        );
    }

    #[test]
    fn giersiepen_increases_with_shear() {
        assert!(
            giersiepen_hi(shear(200.0), duration(1.0)) > giersiepen_hi(shear(100.0), duration(1.0))
        );
    }

    #[test]
    fn giersiepen_increases_with_duration() {
        assert!(
            giersiepen_hi(shear(100.0), duration(2.0)) > giersiepen_hi(shear(100.0), duration(1.0))
        );
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
        let exposure = HemolysisExposure::shear_only(shear(80.0), duration(0.5));
        let expected = giersiepen_hi(shear(80.0), duration(0.5));
        assert!((exposure.compute_index() - expected).abs() < 1e-15);
    }

    #[test]
    fn exposure_with_cavitation_exceeds_base() {
        let base = HemolysisExposure::shear_only(shear(80.0), duration(0.5)).compute_index();
        let cav = HemolysisExposure::new(shear(80.0), duration(0.5), 0.5).compute_index();
        assert!(cav > base);
    }

    // ── taskin_hi ────────────────────────────────────────────────────────────

    #[test]
    fn test_taskin_zero_inputs() {
        assert_eq!(taskin_hi(shear(0.0), duration(1.0)), 0.0);
        assert_eq!(taskin_hi(shear(100.0), duration(0.0)), 0.0);
        assert_eq!(taskin_hi(shear(-10.0), duration(1.0)), 0.0);
        assert_eq!(taskin_hi(shear(10.0), duration(-1.0)), 0.0);
    }

    #[test]
    fn test_taskin_reference_value() {
        // HI = C_T × τ^β_T × t = 1.228e-5 × 100^1.9918 × 1.0
        let expected = TASKIN_C * 100.0_f64.powf(TASKIN_BETA) * 1.0;
        let hi = taskin_hi(shear(100.0), duration(1.0));
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
        assert!(taskin_hi(shear(200.0), duration(1.0)) > taskin_hi(shear(100.0), duration(1.0)));
        // HI increases with exposure time
        assert!(taskin_hi(shear(100.0), duration(2.0)) > taskin_hi(shear(100.0), duration(1.0)));
    }

    #[test]
    fn test_taskin_vs_giersiepen() {
        // At the same inputs, Taskin and Giersiepen should give different values
        // because they use different constants and time dependence.
        let tau = 100.0;
        let t = 1.0;
        let hi_g = giersiepen_hi(shear(tau), duration(t));
        let hi_t = taskin_hi(shear(tau), duration(t));
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
        assert!(v2 > v1, "Higher pressure should give higher jet velocity");
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
        assert_eq!(
            cavitation_hemolysis_amplification(0.0, 1e-6, 101_325.0),
            1.0
        );
        assert_eq!(
            cavitation_hemolysis_amplification(1e-6, 0.0, 101_325.0),
            1.0
        );
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
