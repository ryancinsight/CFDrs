//! Cavitation regime classification: stable vs. inertial cavitation.
//!
//! This module provides tools to classify cavitation behavior based on
//! bubble dynamics, acoustic characteristics, and flow conditions.
//!
//! ## Cavitation Regimes
//!
//! ### Stable Cavitation
//! - Bubbles oscillate about equilibrium radius without collapse
//! - Occurs at lower acoustic pressures or moderate flow velocities
//! - Characterized by sustained oscillations and microstreaming
//! - Less damaging to materials and biological cells
//! - Can produce sonoluminescence under specific conditions
//!
//! ### Inertial Cavitation  
//! - Bubbles grow rapidly and collapse violently
//! - Occurs at high acoustic pressures or high flow velocities
//! - Characterized by transient cavitation events
//! - Produces shock waves, microjets, and high local temperatures
//! - Causes severe material erosion and hemolysis
//!
//! ## Mathematical Criteria
//!
//! ### Blake Threshold
//! ```math
//! P_Blake = P_v + 2σ/R_c * (1 + 2σ/(3R_c(P_∞ - P_v)))
//! ```
//!
//! ### Inertial Cavitation Threshold (Apfel & Holland 1991)
//! ```math
//! P_threshold = P_v + sqrt(8σ/(3R_0)) * (P_∞ + 2σ/R_0)^(1/2)
//! ```

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

use crate::error::Result;

use super::rayleigh_plesset::RayleighPlesset;

/// Cavitation regime types
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum CavitationRegime {
    /// No cavitation occurring
    None,
    /// Stable oscillating cavitation
    Stable,
    /// Transient inertial cavitation
    Inertial,
    /// Mixed regime (both stable and inertial present)
    Mixed,
}

/// Cavitation regime classifier
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationRegimeClassifier<T: RealField + Copy> {
    /// Rayleigh-Plesset bubble model
    bubble_model: RayleighPlesset<T>,
    /// Ambient pressure (Pa)
    ambient_pressure: T,
    /// Acoustic pressure amplitude (Pa), if applicable
    acoustic_pressure: Option<T>,
    /// Acoustic frequency (Hz), if applicable
    acoustic_frequency: Option<T>,
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
        // Check Blake threshold for cavitation inception
        let blake_threshold = self.blake_threshold();
        let current_pressure = if let Some(p_ac) = self.acoustic_pressure {
            self.ambient_pressure - p_ac
        } else {
            self.ambient_pressure
        };

        if current_pressure >= blake_threshold {
            // No cavitation - pressure too high
            return CavitationRegime::None;
        }

        // Calculate inertial threshold
        let inertial_threshold = self.inertial_threshold();

        // Classify based on pressure amplitude
        if let Some(p_ac) = self.acoustic_pressure {
            if p_ac > inertial_threshold {
                CavitationRegime::Inertial
            } else {
                CavitationRegime::Stable
            }
        } else {
            // For hydrodynamic cavitation, use cavitation number
            let sigma = self.cavitation_number(current_pressure);
            
            // Empirical thresholds from literature:
            // σ > 1.5: no cavitation
            // 0.5 < σ < 1.5: stable cavitation
            // σ < 0.5: inertial cavitation
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
        // P_Blake = P_v + (4σ/3R_c)
        let four_thirds = T::from_f64(4.0 / 3.0).unwrap_or_else(|| T::one());
        let r_critical = self.bubble_model.blake_critical_radius(self.ambient_pressure);
        
        self.bubble_model.vapor_pressure
            + four_thirds * self.bubble_model.surface_tension / r_critical
    }

    /// Calculate inertial cavitation threshold (Apfel & Holland 1991)
    pub fn inertial_threshold(&self) -> T {
        // P_threshold = sqrt(8σ/(3R_0) * (P_∞ + 2σ/R_0))
        let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one());
        let eight = T::from_f64(8.0).unwrap_or_else(|| T::one());
        
        let r0 = self.bubble_model.initial_radius;
        let sigma = self.bubble_model.surface_tension;
        let p_inf = self.ambient_pressure;
        
        let term1 = eight * sigma / (three * r0);
        let term2 = p_inf + two * sigma / r0;
        
        (term1 * term2).sqrt()
    }

    /// Calculate cavitation number
    pub fn cavitation_number(&self, local_pressure: T) -> T {
        // σ = (P - P_v) / (0.5 * ρ * U²)
        // Simplified version using pressure difference only
        let pressure_diff = local_pressure - self.bubble_model.vapor_pressure;
        let reference_pressure = self.ambient_pressure - self.bubble_model.vapor_pressure;
        
        if reference_pressure > T::zero() {
            pressure_diff / reference_pressure
        } else {
            T::zero()
        }
    }

    /// Estimate maximum bubble radius for given regime
    pub fn estimate_max_radius(&self, regime: CavitationRegime) -> T {
        match regime {
            CavitationRegime::None => self.bubble_model.initial_radius,
            CavitationRegime::Stable => {
                // Stable cavitation: R_max ≈ 1.5 to 3 × R_0
                let two = T::from_f64(2.0).unwrap_or_else(|| T::one());
                two * self.bubble_model.initial_radius
            }
            CavitationRegime::Inertial => {
                // Inertial cavitation: R_max >> R_0 (up to 10-100×)
                let fifty = T::from_f64(50.0).unwrap_or_else(|| T::one());
                fifty * self.bubble_model.initial_radius
            }
            CavitationRegime::Mixed => {
                let ten = T::from_f64(10.0).unwrap_or_else(|| T::one());
                ten * self.bubble_model.initial_radius
            }
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
        
        // MI = P_ac / sqrt(f_MHz)
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
            CavitationRegime::Stable => {
                // Stable cavitation can produce sonoluminescence
                // under specific resonance conditions (SBSL)
                T::from_f64(0.7).unwrap_or_else(|| T::one())
            }
            CavitationRegime::Inertial => {
                // Inertial cavitation always produces sonoluminescence
                // during violent collapse (MBSL - multi-bubble)
                T::from_f64(0.95).unwrap_or_else(|| T::one())
            }
            CavitationRegime::Mixed => {
                T::from_f64(0.80).unwrap_or_else(|| T::one())
            }
        }
    }

    /// Estimate damage potential (0-1 scale)
    pub fn damage_potential(&self) -> T {
        let regime = self.classify_regime();
        
        match regime {
            CavitationRegime::None => T::zero(),
            CavitationRegime::Stable => {
                // Stable cavitation: minimal damage from microstreaming
                T::from_f64(0.1).unwrap_or_else(|| T::zero())
            }
            CavitationRegime::Inertial => {
                // Inertial cavitation: severe damage from collapse
                T::from_f64(1.0).unwrap_or_else(|| T::one())
            }
            CavitationRegime::Mixed => {
                T::from_f64(0.6).unwrap_or_else(|| T::one())
            }
        }
    }

    /// Estimate hemolysis risk for blood flow
    pub fn hemolysis_risk(&self) -> T {
        let regime = self.classify_regime();
        let damage = self.damage_potential();
        
        match regime {
            CavitationRegime::None => T::zero(),
            CavitationRegime::Stable => {
                // Stable cavitation: low hemolysis from microstreaming
                T::from_f64(0.05).unwrap_or_else(|| T::zero())
            }
            CavitationRegime::Inertial => {
                // Inertial cavitation: high hemolysis from shock waves
                // Scale with pressure amplitude if available
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
            CavitationRegime::Mixed => {
                T::from_f64(0.4).unwrap_or_else(|| T::one())
            }
        }
    }
}

/// Cavitation regime analysis results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CavitationRegimeAnalysis<T: RealField + Copy> {
    /// Identified regime
    pub regime: CavitationRegime,
    /// Blake threshold pressure (Pa)
    pub blake_threshold: T,
    /// Inertial threshold pressure amplitude (Pa)
    pub inertial_threshold: T,
    /// Cavitation number
    pub cavitation_number: T,
    /// Mechanical index (for acoustic cavitation)
    pub mechanical_index: T,
    /// Estimated maximum bubble radius (m)
    pub max_bubble_radius: T,
    /// Sonoluminescence probability (0-1)
    pub sonoluminescence_probability: T,
    /// Material damage potential (0-1)
    pub damage_potential: T,
    /// Hemolysis risk (0-1)
    pub hemolysis_risk: T,
}

impl<T: RealField + Copy + FromPrimitive> CavitationRegimeClassifier<T> {
    /// Perform comprehensive regime analysis
    pub fn analyze(&self) -> Result<CavitationRegimeAnalysis<T>> {
        let regime = self.classify_regime();
        let blake_threshold = self.blake_threshold();
        let inertial_threshold = self.inertial_threshold();
        
        let current_pressure = if let Some(p_ac) = self.acoustic_pressure {
            self.ambient_pressure - p_ac
        } else {
            self.ambient_pressure
        };
        
        let cavitation_number = self.cavitation_number(current_pressure);
        let mechanical_index = self.mechanical_index()?;
        let max_bubble_radius = self.estimate_max_radius(regime);
        let sonoluminescence_probability = self.sonoluminescence_probability();
        let damage_potential = self.damage_potential();
        let hemolysis_risk = self.hemolysis_risk();
        
        Ok(CavitationRegimeAnalysis {
            regime,
            blake_threshold,
            inertial_threshold,
            cavitation_number,
            mechanical_index,
            max_bubble_radius,
            sonoluminescence_probability,
            damage_potential,
            hemolysis_risk,
        })
    }
}

impl<T: RealField + Copy + std::fmt::Display> std::fmt::Display for CavitationRegimeAnalysis<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Cavitation Regime Analysis:\n\
             Regime: {:?}\n\
             Blake Threshold: {} Pa\n\
             Inertial Threshold: {} Pa\n\
             Cavitation Number: {}\n\
             Mechanical Index: {}\n\
             Max Bubble Radius: {} μm\n\
             Sonoluminescence Probability: {:. 1}%\n\
             Damage Potential: {:.1}%\n\
             Hemolysis Risk: {:.1}%",
            self.regime,
            self.blake_threshold,
            self.inertial_threshold,
            self.cavitation_number,
            self.mechanical_index,
            self.max_bubble_radius,
            self.sonoluminescence_probability,
            self.damage_potential,
            self.hemolysis_risk,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn create_test_bubble() -> RayleighPlesset<f64> {
        RayleighPlesset {
            initial_radius: 1e-6,
            liquid_density: 998.0,
            liquid_viscosity: 1.002e-3,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        }
    }

    #[test]
    fn test_regime_classification_none() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            1e5, // High ambient pressure
            None,
            None,
        );
        
        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::None);
    }

    #[test]
    fn test_regime_classification_stable() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            1e5,
            Some(1e4), // Moderate acoustic pressure
            Some(20e3), // 20 kHz
        );
        
        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Stable);
    }

    #[test]
    fn test_regime_classification_inertial() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            1e5,
            Some(1e6), // High acoustic pressure
            Some(20e3),
        );
        
        let regime = classifier.classify_regime();
        assert_eq!(regime, CavitationRegime::Inertial);
    }

    #[test]
    fn test_mechanical_index() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(
            bubble,
            1e5,
            Some(1e6),
            Some(1e6), // 1 MHz
        );
        
        let mi = classifier.mechanical_index().unwrap();
        assert!(mi > 0.0);
        // MI = 1e6 / sqrt(1) = 1e6
        assert!((mi - 1e6).abs() < 1.0);
    }

    #[test]
    fn test_damage_potential_ordering() {
        let bubble = create_test_bubble();
        
        let classifier_none = CavitationRegimeClassifier::new(bubble, 1e5, None, None);
        let classifier_stable = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e4), Some(20e3));
        let classifier_inertial = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(20e3));
        
        let damage_none = classifier_none.damage_potential();
        let damage_stable = classifier_stable.damage_potential();
        let damage_inertial = classifier_inertial.damage_potential();
        
        assert!(damage_none < damage_stable);
        assert!(damage_stable < damage_inertial);
    }

    #[test]
    fn test_hemolysis_risk_inertial_highest() {
        let bubble = create_test_bubble();
        
        let classifier_stable = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e4), Some(20e3));
        let classifier_inertial = CavitationRegimeClassifier::new(bubble, 1e5, Some(1e6), Some(20e3));
        
        let risk_stable = classifier_stable.hemolysis_risk();
        let risk_inertial = classifier_inertial.hemolysis_risk();
        
        assert!(risk_inertial > risk_stable);
        assert!(risk_inertial > 0.5); // Inertial should be high risk
    }

    #[test]
    fn test_full_analysis() {
        let bubble = create_test_bubble();
        let classifier = CavitationRegimeClassifier::new(bubble, 1e5, Some(5e5), Some(20e3));
        
        let analysis = classifier.analyze().unwrap();
        
        assert_eq!(analysis.regime, CavitationRegime::Inertial);
        assert!(analysis.blake_threshold > 0.0);
        assert!(analysis.inertial_threshold > 0.0);
        assert!(analysis.max_bubble_radius > bubble.initial_radius);
        assert!(analysis.sonoluminescence_probability > 0.5);
    }
}
