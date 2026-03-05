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
        Self { shear_pa, duration_s, cavitation_potential }
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
        let expected = GIERSIEPEN_C * 1.0_f64.powf(GIERSIEPEN_ALPHA) * 100.0_f64.powf(GIERSIEPEN_BETA);
        let hi = giersiepen_hi(100.0, 1.0);
        assert!((hi - expected).abs() < 1e-15, "got {hi}, expected {expected}");
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
        assert!((cavitation_amplified_hi(base, 2.0) - cavitation_amplified_hi(base, 1.0)).abs() < 1e-15);
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
}
