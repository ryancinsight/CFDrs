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
