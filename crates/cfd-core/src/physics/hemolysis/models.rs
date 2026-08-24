//! Hemolysis model definitions and damage index calculations.
//!
//! # Theorem — Giersiepen-Wurzinger Power Law (Giersiepen et al. 1990)
//!
//! The fractional increase in plasma-free haemoglobin (Hemolysis Index) for
//! red blood cells exposed to a scalar shear stress $\tau$ for duration $t$
//! follows the empirical power law:
//!
//! ```text
//! HI = ΔHb / Hb = C · τ^α · t^β
//! ```
//!
//! with standard constants $C = 3.62 \times 10^{-7}$, $\alpha = 2.416$,
//! $\beta = 0.785$ (laminar Couette-flow calibration).
//!
//! **Proof sketch.** Giersiepen et al. performed controlled Couette-viscometer
//! experiments varying shear stress (40–255 Pa) and exposure time (7–700 ms).
//! Least-squares regression of log(HI) against log(τ) and log(t) yields the
//! power-law exponents. The model assumes cumulative sub-lethal damage
//! (Leverett et al. 1972) and is valid for $\tau < 800$ Pa and $t < 1$ s.
//!
//! ## References
//!
//! - Giersiepen, M. et al. (1990). "Estimation of shear stress-related blood
//!   damage in heart valve prostheses." *Int. J. Artif. Organs* 13:300–306.
//! - Leverett, L. B. et al. (1972). "Red blood cell damage by shear stress."
//!   *Biophys. J.* 12:257–273.
//! - Heuser, G. & Opitz, R. (1980). "A Couette viscometer for short time
//!   shearing of blood." *Biorheology* 17:17–24.

use crate::error::{Error, Result};
use serde::{Deserialize, Serialize};

// ── Named constants for the Giersiepen power law ─────────────────────────────

/// Giersiepen (1990) calibration constant C for the haemolysis index
/// `HI = C · τ^α · t^β` (τ in Pa, t in s).
///
/// **Reference**: Giersiepen M. et al. (1990). "Estimation of shear stress-related
/// blood damage in heart valve prostheses." *Int. J. Artif. Organs* 13(5):300–306.
pub const GIERSIEPEN_MILLIFLUIDIC_C: f64 = 3.62e-7;

/// Giersiepen (1990) shear-stress exponent α.
///
/// Calibrated from Couette-viscometer experiments at τ ∈ [40, 255] Pa.
pub const GIERSIEPEN_MILLIFLUIDIC_STRESS: f64 = 2.416;

/// Giersiepen (1990) time exponent β.
///
/// Calibrated from Couette-viscometer experiments at t ∈ [7, 700] ms.
pub const GIERSIEPEN_MILLIFLUIDIC_TIME: f64 = 0.785;

/// Conservative cavitation amplification slope for SDT millifluidic devices.
///
/// `HI_amp = base_hi · (1 + CAVITATION_HI_SLOPE · σ.clamp(0, 1))`.
/// At full cavitation potential (σ = 1), this yields 4× the baseline HI.
/// Acoustic bubble collapse generates micro-jets and shockwaves that damage
/// RBC membranes independently of macroscopic steady shear stress.
/// The 3× slope is conservative; experimental values range 2–5×
/// (Brujan 2011 §7.4) depending on bubble proximity and collapse geometry.
pub const CAVITATION_HI_SLOPE: f64 = 3.0;

/// Hemolysis model types
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum HemolysisModel {
    /// Power law model (Giersiepen et al. 1990)
    /// Most widely used for cardiovascular devices
    PowerLaw {
        /// Coefficient C (dimensionless)
        coefficient: f64,
        /// Shear stress exponent α (typically 0.765-2.416)
        stress_exponent: f64,
        /// Time exponent β (typically 0.547-0.785)
        time_exponent: f64,
    },
    /// Zhang model (2011) for Couette flow
    Zhang {
        /// Material coefficient
        coefficient: f64,
        /// Shear rate exponent
        rate_exponent: f64,
    },
    /// Generic linear-threshold damage model.
    LinearThreshold {
        /// Critical shear stress threshold (Pa)
        threshold: f64,
        /// Damage rate above threshold
        damage_rate: f64,
    },
}

impl Default for HemolysisModel {
    fn default() -> Self {
        Self::PowerLaw {
            coefficient: 3.62e-7,
            stress_exponent: 2.416,
            time_exponent: 0.785,
        }
    }
}

impl HemolysisModel {
    /// Create Giersiepen model with standard constants
    pub fn giersiepen_standard() -> Self {
        Self::PowerLaw {
            coefficient: 3.62e-7,
            stress_exponent: 2.416,
            time_exponent: 0.785,
        }
    }

    /// Heuser & Opitz (1980) Couette-viscometer power law:
    /// `HI = 1.8×10⁻⁶ · τ^2.09 · t^0.765`.
    ///
    /// **Reference**: Heuser, G. & Opitz, R. (1980). "A Couette viscometer for
    /// short time shearing of blood." *Biorheology* 17:17–24.
    pub fn heuser_opitz_couette() -> Self {
        Self::PowerLaw {
            coefficient: 1.8e-6,
            stress_exponent: 2.09,
            time_exponent: 0.765,
        }
    }

    /// Create a generic linear-threshold damage model.
    ///
    /// Note: this is *not* the published Heuser & Opitz (1980) correlation —
    /// that is the power law returned by [`Self::heuser_opitz_couette`].
    pub fn linear_threshold(threshold: f64, damage_rate: f64) -> Self {
        Self::LinearThreshold {
            threshold,
            damage_rate,
        }
    }

    /// Zhang model for Couette flow (Zhang et al. 2011):
    /// `HI = 1.86×10⁻⁴ · γ̇^1.84 · t` with γ̇ the shear rate in s⁻¹.
    pub fn zhang() -> Self {
        Self::Zhang {
            coefficient: 1.86e-4,
            rate_exponent: 1.84,
        }
    }

    /// Giersiepen (1990) model applied to millifluidic and blood-processing
    /// device screening (SDT, micro-pump, oxygenator, venturi).
    ///
    /// Uses the standard Couette calibration — no device-specific
    /// recalibration exists in the literature. Prior to 2026-08-20 this
    /// variant carried fabricated constants (C = 3.62×10⁻⁵, β = 1.991,
    /// α = 0.765) presented as "Giersiepen Table 1"; the actual published
    /// values are C = 3.62×10⁻⁷, α = 2.416, β = 0.785.
    ///
    /// # Theorem — Giersiepen Power-Law (Giersiepen et al. 1990)
    ///
    /// The fractional haemoglobin release (Haemolysis Index) for red blood cells
    /// exposed to scalar shear stress τ for duration t follows the empirical law:
    ///
    /// ```text
    /// HI = C · τ^α · t^β
    /// ```
    ///
    /// with C = 3.62 × 10⁻⁷,  α = 2.416 (shear exponent),  β = 0.785 (time exponent).
    /// These constants were calibrated against Couette-viscometer experiments at
    /// τ ∈ [40, 255] Pa and t ∈ [7, 700] ms (Table 1 of Giersiepen et al. 1990).
    ///
    /// **Proof sketch.** Least-squares regression of log(HI) against log(τ) and
    /// log(t) from controlled Couette-viscometer experiments yields the power-law
    /// exponents.  The model assumes cumulative sub-lethal damage (Leverett 1972)
    /// and is valid for τ < 800 Pa and t < 1 s.
    ///
    /// **References.**
    /// - Giersiepen, M. et al. (1990). "Estimation of shear stress-related blood
    ///   damage in heart valve prostheses." *Int. J. Artif. Organs* 13(5):300–306.
    /// - Leverett, L.B. et al. (1972). "Red blood cell damage by shear stress."
    ///   *Biophys. J.* 12:257–273.
    pub fn giersiepen_millifluidic() -> Self {
        Self::PowerLaw {
            coefficient: GIERSIEPEN_MILLIFLUIDIC_C,
            stress_exponent: GIERSIEPEN_MILLIFLUIDIC_STRESS,
            time_exponent: GIERSIEPEN_MILLIFLUIDIC_TIME,
        }
    }

    /// Amplify a baseline Haemolysis Index by local cavitation potential σ ∈ [0, 1].
    ///
    /// # Theorem — Conservative Cavitation Amplification
    ///
    /// Acoustic bubble collapse (inertial cavitation) generates micro-jets with
    /// localised shear amplification and pressure shockwaves that rupture RBC
    /// membranes independently of macroscopic steady shear.  The amplified HI is:
    ///
    /// ```text
    /// HI_amp = base_hi · (1 + CAVITATION_HI_SLOPE · σ.clamp(0, 1))
    /// ```
    ///
    /// At σ = 1 (full cavitation): HI_amp = 4 · base_hi (4× amplification).
    /// At σ = 0 (no cavitation):   HI_amp = base_hi      (no change).
    ///
    /// **References.**
    /// - Brujan, E.A. (2011). *Cavitation in Non-Newtonian Fluids.* Springer §7.4.
    #[must_use]
    pub fn cavitation_amplified(base_hi: f64, cav_potential: f64) -> f64 {
        base_hi * (1.0 + CAVITATION_HI_SLOPE * cav_potential.clamp(0.0, 1.0))
    }

    /// Calculate blood damage index from shear stress and exposure time
    pub fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> Result<f64> {
        if shear_stress < 0.0 {
            return Err(Error::InvalidConfiguration(
                "Shear stress must be non-negative".to_string(),
            ));
        }
        if exposure_time < 0.0 {
            return Err(Error::InvalidConfiguration(
                "Exposure time must be non-negative".to_string(),
            ));
        }

        match self {
            Self::PowerLaw {
                coefficient,
                stress_exponent,
                time_exponent,
            } => Ok(coefficient
                * shear_stress.powf(*stress_exponent)
                * exposure_time.powf(*time_exponent)),

            Self::Zhang {
                coefficient,
                rate_exponent,
            } => {
                let shear_rate = shear_stress / crate::physics::fluid::blood::constants::INFINITE_SHEAR_VISCOSITY;
                Ok(coefficient * shear_rate.powf(*rate_exponent) * exposure_time)
            }

            Self::LinearThreshold {
                threshold,
                damage_rate,
            } => {
                if shear_stress > *threshold {
                    Ok(damage_rate * (shear_stress - threshold) * exposure_time)
                } else {
                    Ok(0.0)
                }
            }
        }
    }
}
