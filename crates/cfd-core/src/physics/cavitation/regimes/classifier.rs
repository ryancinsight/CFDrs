//! Cavitation regime classifier: stable vs. inertial cavitation.
//!
//! Classifies cavitation behavior based on bubble dynamics, acoustic
//! characteristics, and flow conditions using Blake and Apfel–Holland
//! thresholds.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

use crate::error::Result;

use super::super::rayleigh_plesset::RayleighPlesset;
use super::types::CavitationRegime;

/// Cavitation regime classifier
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationRegimeClassifier<T: RealField + Copy> {
    /// Rayleigh-Plesset bubble model
    pub bubble_model: RayleighPlesset<T>,
    /// Ambient pressure (Pa)
    pub ambient_pressure: T,
    /// Acoustic pressure amplitude (Pa), if applicable
    pub acoustic_pressure: Option<T>,
    /// Acoustic frequency (Hz), if applicable
    pub acoustic_frequency: Option<T>,
}

impl<T: RealField + Copy + FromPrimitive> CavitationRegimeClassifier<T> {
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
            let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
            let v_sq =
                (self.ambient_pressure - min_pressure) * two / self.bubble_model.liquid_density;
            let current_velocity = if v_sq > T::zero() {
                v_sq.sqrt()
            } else {
                T::from_f64(1e-6).unwrap_or_else(|| T::one())
            };

            let sigma = self.cavitation_number(min_pressure, current_velocity);

            let one_half = T::from_f64(0.5).unwrap_or_else(|| T::one());
            let one_and_half = T::from_f64(1.5).unwrap_or_else(|| T::one());

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
        let four_thirds = T::from_f64(4.0 / 3.0).unwrap_or_else(|| T::one());
        let r_critical = self
            .bubble_model
            .blake_critical_radius(self.ambient_pressure);

        self.bubble_model.vapor_pressure
            + four_thirds * self.bubble_model.surface_tension / r_critical
    }

    /// Calculate inertial cavitation threshold (Apfel & Holland 1991)
    pub fn inertial_threshold(&self) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one());
        let eight = T::from_f64(8.0).unwrap_or_else(|| T::one());

        let r0 = self.bubble_model.initial_radius;
        let sigma = self.bubble_model.surface_tension;
        let p_inf = self.ambient_pressure;

        let inner = (eight * sigma) / (three * r0) * (p_inf + two * sigma / r0);
        if inner > T::zero() {
            inner.sqrt()
        } else {
            T::zero()
        }
    }

    /// Estimate maximum bubble radius for a given cavitation regime
    pub(super) fn estimate_max_radius(&self, regime: CavitationRegime) -> T {
        match regime {
            CavitationRegime::None => self.bubble_model.initial_radius,
            CavitationRegime::Stable => {
                let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
                two * self.bubble_model.initial_radius
            }
            CavitationRegime::Inertial => {
                let fifty = T::from_f64(50.0).unwrap_or_else(|| T::one());
                fifty * self.bubble_model.initial_radius
            }
            CavitationRegime::Mixed => {
                let ten = T::from_f64(10.0).unwrap_or_else(|| T::one());
                ten * self.bubble_model.initial_radius
            }
        }
    }

    /// Calculate cavitation number σ = (P - P_v) / (0.5 ρ v²)
    pub(super) fn cavitation_number(&self, local_pressure: T, velocity: T) -> T {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one());
        let dynamic_pressure = half * self.bubble_model.liquid_density * velocity * velocity;
        if dynamic_pressure > T::default_epsilon() {
            (local_pressure - self.bubble_model.vapor_pressure) / dynamic_pressure
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one())
        }
    }

    /// Calculate mechanical index (MI) for ultrasound cavitation
    pub fn mechanical_index(&self) -> Result<T> {
        let Some(p_ac) = self.acoustic_pressure else {
            return Ok(T::zero());
        };
        let Some(freq) = self.acoustic_frequency else {
            return Ok(T::zero());
        };

        let mhz = T::from_f64(1e6).unwrap_or_else(|| T::one());
        let freq_mhz = freq / mhz;

        if freq_mhz > T::zero() {
            Ok(p_ac / freq_mhz.sqrt())
        } else {
            Ok(T::zero())
        }
    }

    /// Estimate sonoluminescence probability based on regime
    pub fn sonoluminescence_probability(&self) -> T {
        let regime = self.classify_regime();

        match regime {
            CavitationRegime::None => T::zero(),
            CavitationRegime::Stable => T::from_f64(0.7).unwrap_or_else(|| T::one()),
            CavitationRegime::Inertial => T::from_f64(0.95).unwrap_or_else(|| T::one()),
            CavitationRegime::Mixed => T::from_f64(0.80).unwrap_or_else(|| T::one()),
        }
    }

    /// Estimate damage potential (0-1 scale)
    pub fn damage_potential(&self) -> T {
        let regime = self.classify_regime();

        match regime {
            CavitationRegime::None => T::zero(),
            CavitationRegime::Stable => T::from_f64(0.1).unwrap_or_else(|| T::zero()),
            CavitationRegime::Inertial => T::from_f64(1.0).unwrap_or_else(|| T::one()),
            CavitationRegime::Mixed => T::from_f64(0.6).unwrap_or_else(|| T::one()),
        }
    }

    /// Estimate hemolysis risk for blood flow
    pub fn hemolysis_risk(&self) -> T {
        let regime = self.classify_regime();
        let damage = self.damage_potential();

        match regime {
            CavitationRegime::None => T::zero(),
            CavitationRegime::Stable => T::from_f64(0.05).unwrap_or_else(|| T::zero()),
            CavitationRegime::Inertial => {
                if let Some(p_ac) = self.acoustic_pressure {
                    let threshold = self.inertial_threshold();
                    if threshold > T::zero() {
                        (p_ac / threshold).min(T::one())
                    } else {
                        damage
                    }
                } else {
                    damage
                }
            }
            CavitationRegime::Mixed => T::from_f64(0.4).unwrap_or_else(|| T::one()),
        }
    }
}
