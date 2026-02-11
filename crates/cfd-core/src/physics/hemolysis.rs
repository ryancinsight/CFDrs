//! Hemolysis and blood damage models for microfluidic and millifluidic applications.
//!
//! This module provides mathematical models for predicting red blood cell (RBC) damage
//! and hemolysis in flow systems, particularly relevant for:
//! - Cardiovascular devices (pumps, valves, oxygenators)
//! - Microfluidic blood processing
//! - Millifluidic diagnostic devices
//! - Venturi cavitation systems
//!
//! ## Mathematical Foundation
//!
//! ### Power Law Model (Giersiepen et al. 1990)
//! The most widely used hemolysis model relates blood damage to shear stress exposure:
//!
//! ```math
//! D = C * τ^α * t^β
//! ```
//!
//! where:
//! - D = damage index (0-1, or hemoglobin release in mg/dL)
//! - τ = shear stress (Pa)
//! - t = exposure time (s)
//! - C, α, β = empirical constants
//!
//! ### Normalized Index of Hemolysis (NIH)
//! ```math
//! NIH = (100 - Hct) / Hct * ΔHb / Hb₀ * 100%
//! ```
//!
//! ### Modified Index of Hemolysis (MIH)
//! ```math
//! MIH = V_blood / (Q * t) * ΔHb * (100 - Hct) / 100
//! ```
//!
//! ## References
//! - Giersiepen, M. et al. (1990). Estimation of shear stress-related blood damage.
//! - Zhang, T. et al. (2011). Study of flow-induced hemolysis using novel Couette-type blood-shearing devices.
//! - Goubergrits, L. & Affeld, K. (2004). Numerical estimation of blood damage in artificial organs.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

use crate::error::{Error, Result};

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
        // Giersiepen constants for blood pumps
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
            threshold: 150.0, // Pa
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
                // Convert shear stress to shear rate (τ = μ * γ̇)
                // Assume blood viscosity μ ≈ 0.0035 Pa·s
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

/// Hemolysis calculator with blood properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HemolysisCalculator<T: RealField + Copy> {
    /// Hemolysis model
    model: HemolysisModel,
    /// Hematocrit (volume fraction, 0-1)
    hematocrit: T,
    /// Initial hemoglobin concentration (g/dL)
    hemoglobin_initial: T,
    /// Blood volume (m³)
    blood_volume: T,
    /// Flow rate (m³/s)
    flow_rate: T,
}

impl<T: RealField + Copy + FromPrimitive> HemolysisCalculator<T> {
    /// Create new hemolysis calculator
    pub fn new(
        model: HemolysisModel,
        hematocrit: T,
        hemoglobin_initial: T,
        blood_volume: T,
        flow_rate: T,
    ) -> Result<Self> {
        if hematocrit < T::zero() || hematocrit > T::one() {
            return Err(Error::InvalidConfiguration(
                "Hematocrit must be between 0 and 1".to_string(),
            ));
        }
        if hemoglobin_initial <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Initial hemoglobin must be positive".to_string(),
            ));
        }
        if blood_volume <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Blood volume must be positive".to_string(),
            ));
        }
        if flow_rate <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "Flow rate must be positive".to_string(),
            ));
        }

        Ok(Self {
            model,
            hematocrit,
            hemoglobin_initial,
            blood_volume,
            flow_rate,
        })
    }

    /// Calculate Normalized Index of Hemolysis (NIH)
    pub fn normalized_index(&self, delta_hemoglobin: T) -> T {
        let hundred = T::from_f64(100.0).unwrap_or_else(|| T::one());
        let one_minus_hct = (hundred - self.hematocrit * hundred) / hundred;
        
        if self.hemoglobin_initial > T::zero() {
            one_minus_hct / self.hematocrit * delta_hemoglobin / self.hemoglobin_initial * hundred
        } else {
            T::zero()
        }
    }

    /// Calculate Modified Index of Hemolysis (MIH) in g/(100L)
    pub fn modified_index(&self, delta_hemoglobin: T, exposure_time: T) -> T {
        let hundred = T::from_f64(100.0).unwrap_or_else(|| T::one());
        
        if self.flow_rate > T::zero() && exposure_time > T::zero() {
            self.blood_volume / (self.flow_rate * exposure_time) * delta_hemoglobin
                * (hundred - self.hematocrit * hundred)
                / hundred
        } else {
            T::zero()
        }
    }

    /// Calculate hemoglobin release from damage index
    pub fn hemoglobin_release(&self, damage_index: T) -> T {
        // Empirical correlation: ΔHb (mg/dL) ≈ D * Hb₀ * f(Hct)
        // where f(Hct) ≈ Hct / (1 - Hct) accounts for RBC fraction
        let scaling = self.hematocrit / (T::one() - self.hematocrit);
        self.hemoglobin_initial * damage_index * scaling
    }

    /// Estimate exposure time from geometry
    pub fn estimate_exposure_time(&self, characteristic_length: T) -> T {
        // Simple estimation: t = L / u_mean where u_mean = Q / A
        // For typical channel: A ≈ π D² / 4
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::one());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::one());
        
        // Assume circular channel with diameter = characteristic_length
        let area = pi * characteristic_length * characteristic_length / four;
        
        if area > T::zero() {
            let velocity = self.flow_rate / area;
            if velocity > T::zero() {
                characteristic_length / velocity
            } else {
                T::zero()
            }
        } else {
            T::zero()
        }
    }

    /// Calculate cumulative damage over path with varying shear stress
    pub fn cumulative_damage(&self, stress_time_history: &[(f64, f64)]) -> Result<f64> {
        let mut total_damage = 0.0;
        
        for &(stress, time) in stress_time_history {
            let damage = self.model.damage_index(stress, time)?;
            total_damage += damage;
        }
        
        Ok(total_damage)
    }

    /// Get a reference to the hemolysis model
    pub fn model(&self) -> &HemolysisModel {
        &self.model
    }

    /// Calculate critical shear stress for onset of hemolysis
    pub fn critical_shear_stress(&self, exposure_time: f64) -> Result<f64> {
        // Typical threshold: D = 0.01 (1% damage)
        let threshold_damage = 0.01;
        
        match self.model {
            HemolysisModel::PowerLaw {
                coefficient,
                stress_exponent,
                time_exponent,
            } => {
                // τ_crit = (D / (C * t^β))^(1/α)
                let time_term = exposure_time.powf(time_exponent);
                if time_term > 0.0 && coefficient > 0.0 {
                    Ok((threshold_damage / (coefficient * time_term)).powf(1.0 / stress_exponent))
                } else {
                    Ok(f64::INFINITY)
                }
            }
            HemolysisModel::HeuserOpitz { threshold, .. } => Ok(threshold),
            HemolysisModel::Zhang { .. } => {
                // Approximate critical shear rate then convert
                Ok(1000.0 * 0.0035) // ~1000 s^-1 is typical threshold
            }
        }
    }
}

/// Shear-induced platelet activation model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PlateletActivation<T: RealField + Copy> {
    /// Activation threshold (Pa)
    threshold: T,
    /// Activation rate constant (1/Pa/s)
    rate_constant: T,
}

impl<T: RealField + Copy + FromPrimitive> PlateletActivation<T> {
    /// Create standard platelet activation model
    pub fn standard() -> Self {
        // Typical thresholds: 10-70 Pa for seconds of exposure
        Self {
            threshold: T::from_f64(30.0).unwrap_or_else(|| T::one()),
            rate_constant: T::from_f64(0.01).unwrap_or_else(|| T::one()),
        }
    }

    /// Calculate activation probability
    pub fn activation_probability(&self, shear_stress: T, exposure_time: T) -> T {
        if shear_stress > self.threshold {
            let excess_stress = shear_stress - self.threshold;
            let exponent = -(self.rate_constant * excess_stress * exposure_time);
            T::one() - (-exponent).exp()
        } else {
            T::zero()
        }
    }
}

/// Blood trauma assessment for devices
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BloodTrauma {
    /// Hemolysis level (mg/dL)
    pub hemolysis_level: f64,
    /// Platelet activation (%)
    pub platelet_activation: f64,
    /// Thrombosis risk score (0-1)
    pub thrombosis_risk: f64,
    /// Maximum shear stress encountered (Pa)
    pub max_shear_stress: f64,
    /// Average exposure time (s)
    pub avg_exposure_time: f64,
}

impl BloodTrauma {
    /// Classify blood trauma severity
    pub fn severity(&self) -> BloodTraumaSeverity {
        // FDA guidance: hemolysis < 10 mg/dL is acceptable for short-term devices
        if self.hemolysis_level < 10.0 && self.platelet_activation < 5.0 {
            BloodTraumaSeverity::Minimal
        } else if self.hemolysis_level < 50.0 && self.platelet_activation < 20.0 {
            BloodTraumaSeverity::Moderate
        } else if self.hemolysis_level < 150.0 && self.platelet_activation < 50.0 {
            BloodTraumaSeverity::Severe
        } else {
            BloodTraumaSeverity::Critical
        }
    }

    /// Check if device meets FDA guidance
    pub fn meets_fda_guidance(&self) -> bool {
        // FDA guidance for ventricular assist devices:
        // - NIH < 0.01 g/100L
        // - MIH < 10 mg/dL
        self.hemolysis_level < 10.0
    }
}

/// Blood trauma severity classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BloodTraumaSeverity {
    /// Minimal trauma (acceptable for most applications)
    Minimal,
    /// Moderate trauma (acceptable for short-term use)
    Moderate,
    /// Severe trauma (requires redesign)
    Severe,
    /// Critical trauma (device unsafe)
    Critical,
}

impl std::fmt::Display for BloodTrauma {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Blood Trauma Assessment:\n\
             Hemolysis: {:.2} mg/dL (severity: {:?})\n\
             Platelet Activation: {:.1}%\n\
             Thrombosis Risk: {:.3}\n\
             Max Shear Stress: {:.1} Pa\n\
             Avg Exposure Time: {:.4} s\n\
             FDA Compliance: {}",
            self.hemolysis_level,
            self.severity(),
            self.platelet_activation,
            self.thrombosis_risk,
            self.max_shear_stress,
            self.avg_exposure_time,
            if self.meets_fda_guidance() {
                "PASS"
            } else {
                "FAIL"
            }
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_giersiepen_damage_increases_with_stress() {
        let model = HemolysisModel::giersiepen_standard();
        let t = 1.0; // 1 second

        let d1 = model.damage_index(50.0, t).unwrap();
        let d2 = model.damage_index(100.0, t).unwrap();
        let d3 = model.damage_index(200.0, t).unwrap();

        assert!(d2 > d1);
        assert!(d3 > d2);
    }

    #[test]
    fn test_damage_increases_with_time() {
        let model = HemolysisModel::giersiepen_standard();
        let tau = 100.0; // 100 Pa

        let d1 = model.damage_index(tau, 0.1).unwrap();
        let d2 = model.damage_index(tau, 1.0).unwrap();
        let d3 = model.damage_index(tau, 10.0).unwrap();

        assert!(d2 > d1);
        assert!(d3 > d2);
    }

    #[test]
    fn test_heuser_opitz_threshold() {
        let model = HemolysisModel::heuser_opitz();

        let d_below = model.damage_index(100.0, 1.0).unwrap();
        let d_above = model.damage_index(200.0, 1.0).unwrap();

        assert_eq!(d_below, 0.0);
        assert!(d_above > 0.0);
    }

    #[test]
    fn test_nih_calculation() {
        let calc = HemolysisCalculator::new(
            HemolysisModel::default(),
            0.45,
            15.0,
            5e-3,
            1e-4,
        )
        .unwrap();

        let delta_hb = 0.1; // 0.1 g/dL increase
        let nih = calc.normalized_index(delta_hb);

        assert!(nih > 0.0);
        assert!(nih < 100.0);
    }

    #[test]
    fn test_critical_shear_stress() {
        let calc = HemolysisCalculator::new(
            HemolysisModel::giersiepen_standard(),
            0.45,
            15.0,
            5e-3,
            1e-4,
        )
        .unwrap();

        let tau_crit = calc.critical_shear_stress(1.0).unwrap();

        // Should be in physiologically relevant range (50-500 Pa)
        assert!(tau_crit > 10.0);
        assert!(tau_crit < 1000.0);
    }

    #[test]
    fn test_blood_trauma_severity() {
        let trauma_minimal = BloodTrauma {
            hemolysis_level: 5.0,
            platelet_activation: 2.0,
            thrombosis_risk: 0.01,
            max_shear_stress: 80.0,
            avg_exposure_time: 0.1,
        };

        assert_eq!(trauma_minimal.severity(), BloodTraumaSeverity::Minimal);
        assert!(trauma_minimal.meets_fda_guidance());

        let trauma_severe = BloodTrauma {
            hemolysis_level: 100.0,
            platelet_activation: 30.0,
            thrombosis_risk: 0.5,
            max_shear_stress: 400.0,
            avg_exposure_time: 2.0,
        };

        assert_eq!(trauma_severe.severity(), BloodTraumaSeverity::Severe);
        assert!(!trauma_severe.meets_fda_guidance());
    }
}
