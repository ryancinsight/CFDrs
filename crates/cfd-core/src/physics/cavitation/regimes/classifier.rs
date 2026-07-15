//! Cavitation regime classifier: stable vs. inertial cavitation.
//!
//! Classifies cavitation behavior based on bubble dynamics, acoustic
//! characteristics, and flow conditions using Blake and Apfel–Holland
//! thresholds.

use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

use crate::error::Result;

use super::super::rayleigh_plesset::RayleighPlesset;
use super::types::CavitationRegime;

/// Cavitation regime classifier
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationRegimeClassifier<T: FloatElement + Copy> {
    /// Rayleigh-Plesset bubble model
    pub(super) bubble_model: RayleighPlesset<T>,
    /// Ambient pressure (Pa)
    pub(super) ambient_pressure: T,
    /// Acoustic pressure amplitude (Pa), if applicable
    pub(super) acoustic_pressure: Option<T>,
    /// Acoustic frequency (Hz), if applicable
    pub(super) acoustic_frequency: Option<T>,
}

impl<T: FloatElement + Copy> CavitationRegimeClassifier<T> {
    /// Create new cavitation regime classifier
    pub fn new(
        bubble_model: RayleighPlesset<T>,
        ambient_pressure: T,
        acoustic_pressure: Option<T>,
        acoustic_frequency: Option<T>,
    ) -> Self {
        Self {
            bubble_model,
            ambient_pressure,
            acoustic_pressure,
            acoustic_frequency,
        }
    }

    /// Classify cavitation regime based on bubble behavior
    pub fn classify_regime(&self) -> CavitationRegime {
        let blake_threshold = self.blake_threshold();

        let min_pressure = if let Some(p_ac) = self.acoustic_pressure {
            self.ambient_pressure - p_ac
        } else {
            self.ambient_pressure
        };

        if min_pressure >= blake_threshold {
            return CavitationRegime::None;
        }

        let inertial_threshold = self.inertial_threshold();

        if let Some(p_ac) = self.acoustic_pressure {
            if p_ac > inertial_threshold {
                CavitationRegime::Inertial
            } else {
                CavitationRegime::Stable
            }
        } else {
            let two = <T as FloatElement>::from_f64(2.0);
            let v_sq =
                (self.ambient_pressure - min_pressure) * two / self.bubble_model.liquid_density;
            let current_velocity = if v_sq > <T as NumericElement>::ZERO {
                <T as NumericElement>::sqrt(v_sq)
            } else {
                <T as FloatElement>::from_f64(1e-6)
            };

            let sigma = self.cavitation_number(min_pressure, current_velocity);

            let one_half = <T as FloatElement>::from_f64(0.5);
            let one_and_half = <T as FloatElement>::from_f64(1.5);

            if sigma > one_and_half {
                CavitationRegime::None
            } else if sigma > one_half {
                CavitationRegime::Stable
            } else {
                CavitationRegime::Inertial
            }
        }
    }

    /// Calculate Blake threshold pressure
    pub fn blake_threshold(&self) -> T {
        let four_thirds = <T as FloatElement>::from_f64(4.0 / 3.0);
        let r_critical = self
            .bubble_model
            .blake_critical_radius(self.ambient_pressure);

        self.bubble_model.vapor_pressure
            + four_thirds * self.bubble_model.surface_tension / r_critical
    }

    /// Calculate inertial cavitation threshold (Apfel & Holland 1991)
    pub fn inertial_threshold(&self) -> T {
        let two = <T as FloatElement>::from_f64(2.0);
        let three = <T as FloatElement>::from_f64(3.0);
        let eight = <T as FloatElement>::from_f64(8.0);

        let r0 = self.bubble_model.initial_radius;
        let sigma = self.bubble_model.surface_tension;
        let p_inf = self.ambient_pressure;

        let inner = (eight * sigma) / (three * r0) * (p_inf + two * sigma / r0);
        if inner > <T as NumericElement>::ZERO {
            <T as NumericElement>::sqrt(inner)
        } else {
            <T as NumericElement>::ZERO
        }
    }

    /// Estimate maximum bubble radius for a given cavitation regime
    pub(super) fn estimate_max_radius(&self, regime: CavitationRegime) -> T {
        match regime {
            CavitationRegime::None => self.bubble_model.initial_radius,
            CavitationRegime::Stable => {
                let two = <T as FloatElement>::from_f64(2.0);
                two * self.bubble_model.initial_radius
            }
            CavitationRegime::Inertial => {
                let fifty = <T as FloatElement>::from_f64(50.0);
                fifty * self.bubble_model.initial_radius
            }
            CavitationRegime::Mixed => {
                let ten = <T as FloatElement>::from_f64(10.0);
                ten * self.bubble_model.initial_radius
            }
        }
    }

    /// Calculate cavitation number σ = (P - P_v) / (0.5 ρ v²)
    pub(super) fn cavitation_number(&self, local_pressure: T, velocity: T) -> T {
        let half = <T as FloatElement>::from_f64(0.5);
        let dynamic_pressure = half * self.bubble_model.liquid_density * velocity * velocity;
        if dynamic_pressure > <T as FloatElement>::from_f64(f64::EPSILON) {
            (local_pressure - self.bubble_model.vapor_pressure) / dynamic_pressure
        } else {
            <T as FloatElement>::from_f64(1e10)
        }
    }

    /// Calculate mechanical index (MI) for ultrasound cavitation
    pub fn mechanical_index(&self) -> Result<T> {
        let Some(p_ac) = self.acoustic_pressure else {
            return Ok(<T as NumericElement>::ZERO);
        };
        let Some(freq) = self.acoustic_frequency else {
            return Ok(<T as NumericElement>::ZERO);
        };

        let mhz = <T as FloatElement>::from_f64(1e6);
        let freq_mhz = freq / mhz;

        if freq_mhz > <T as NumericElement>::ZERO {
            Ok(p_ac / <T as NumericElement>::sqrt(freq_mhz))
        } else {
            Ok(<T as NumericElement>::ZERO)
        }
    }

    /// Estimate sonoluminescence probability based on regime
    pub fn sonoluminescence_probability(&self) -> T {
        let regime = self.classify_regime();

        match regime {
            CavitationRegime::None => <T as NumericElement>::ZERO,
            CavitationRegime::Stable => <T as FloatElement>::from_f64(0.7),
            CavitationRegime::Inertial => <T as FloatElement>::from_f64(0.95),
            CavitationRegime::Mixed => <T as FloatElement>::from_f64(0.80),
        }
    }

    /// Estimate damage potential (0-1 scale)
    pub fn damage_potential(&self) -> T {
        let regime = self.classify_regime();

        match regime {
            CavitationRegime::None => <T as NumericElement>::ZERO,
            CavitationRegime::Stable => <T as FloatElement>::from_f64(0.1),
            CavitationRegime::Inertial => <T as NumericElement>::ONE,
            CavitationRegime::Mixed => <T as FloatElement>::from_f64(0.6),
        }
    }

    /// Estimate hemolysis risk for blood flow
    pub fn hemolysis_risk(&self) -> T {
        let regime = self.classify_regime();
        let damage = self.damage_potential();

        match regime {
            CavitationRegime::None => <T as NumericElement>::ZERO,
            CavitationRegime::Stable => <T as FloatElement>::from_f64(0.05),
            CavitationRegime::Inertial => {
                if let Some(p_ac) = self.acoustic_pressure {
                    let threshold = self.inertial_threshold();
                    if threshold > <T as NumericElement>::ZERO {
                        (p_ac / threshold).min_scalar(<T as NumericElement>::ONE)
                    } else {
                        damage
                    }
                } else {
                    damage
                }
            }
            CavitationRegime::Mixed => <T as FloatElement>::from_f64(0.4),
        }
    }
}
