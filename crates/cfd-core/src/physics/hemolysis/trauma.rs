//! Blood trauma assessment and platelet activation models.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

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
