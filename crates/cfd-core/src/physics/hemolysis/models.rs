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
//! with standard constants $C = 3.62 \times 10^{-5}$, $\alpha = 2.416$,
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

// ── Named constants for millifluidic Giersiepen model ─────────────────────────

/// Giersiepen (1990) Table 1: Couette viscometer calibration constant C for
/// blood processing and millifluidic devices.  `HI = C · τ^β · t^α`.
///
/// **Reference**: Giersiepen M. et al. (1990). "Estimation of shear stress-related
/// blood damage in heart valve prostheses." *Int. J. Artif. Organs* 13(5):300–306.
pub const GIERSIEPEN_MILLIFLUIDIC_C: f64 = 3.62e-5;

/// Giersiepen (1990) shear-stress exponent β for millifluidic devices.
///
/// Calibrated from Couette-viscometer experiments at τ ∈ [40, 255] Pa.
pub const GIERSIEPEN_MILLIFLUIDIC_STRESS: f64 = 1.991;

/// Giersiepen (1990) time exponent α for millifluidic devices.
///
/// Calibrated from Couette-viscometer experiments at t ∈ [7, 700] ms.
pub const GIERSIEPEN_MILLIFLUIDIC_TIME: f64 = 0.765;

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
    /// Heuser-Opitz model (1980)
    HeuserOpitz {
        /// Critical shear stress threshold (Pa)
        threshold: f64,
        /// Damage rate above threshold
        damage_rate: f64,
    },
}

impl Default for HemolysisModel {
    fn default() -> Self {
        Self::PowerLaw {
            coefficient: 3.62e-5,
            stress_exponent: 2.416,
            time_exponent: 0.785,
        }
    }
}

impl HemolysisModel {
    /// Create Giersiepen model with standard constants
    pub fn giersiepen_standard() -> Self {
        Self::PowerLaw {
            coefficient: 3.62e-5,
            stress_exponent: 2.416,
            time_exponent: 0.785,
        }
    }

    /// Create Giersiepen model for turbulent flow
    pub fn giersiepen_turbulent() -> Self {
        Self::PowerLaw {
            coefficient: 1.8e-6,
            stress_exponent: 1.991,
            time_exponent: 0.765,
        }
    }

    /// Create Giersiepen model for laminar flow
    pub fn giersiepen_laminar() -> Self {
        Self::PowerLaw {
            coefficient: 1.228e-5,
            stress_exponent: 1.9918,
            time_exponent: 0.6606,
        }
    }

    /// Create Zhang model for Couette flow
    pub fn zhang() -> Self {
        Self::Zhang {
            coefficient: 1.86e-4,
            rate_exponent: 1.84,
        }
    }

    /// Create Heuser-Opitz threshold model
    pub fn heuser_opitz() -> Self {
        Self::HeuserOpitz {
            threshold: 150.0,
            damage_rate: 0.01,
        }
    }

    /// Giersiepen (1990) model validated for millifluidic and blood-processing
    /// devices (SDT, micro-pump, oxygenator, venturi).
    ///
    /// # Theorem — Giersiepen Power-Law (Giersiepen et al. 1990)
    ///
    /// The fractional haemoglobin release (Haemolysis Index) for red blood cells
    /// exposed to scalar shear stress τ for duration t follows the empirical law:
    ///
    /// ```text
    /// HI = C · τ^β · t^α
    /// ```
    ///
    /// with C = 3.62 × 10⁻⁵,  β = 1.991 (shear exponent),  α = 0.765 (time exponent).
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
            coefficient:     GIERSIEPEN_MILLIFLUIDIC_C,
            stress_exponent: GIERSIEPEN_MILLIFLUIDIC_STRESS,
            time_exponent:   GIERSIEPEN_MILLIFLUIDIC_TIME,
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
                let shear_rate = shear_stress / 0.0035;
                Ok(coefficient * shear_rate.powf(*rate_exponent) * exposure_time)
            }

            Self::HeuserOpitz {
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
