use super::{SerpentineGeometry, SerpentineMixingSolution};
use crate::scalar::from_f64;
use crate::scalar::Cfd2dScalar;
use eunomia::FloatElement;
use serde::{Deserialize, Serialize};

/// Validator for serpentine mixing channels
pub struct SerpentineValidator<T: Cfd2dScalar + Copy> {
    geometry: SerpentineGeometry<T>,
}

impl<T: Cfd2dScalar + Copy + FloatElement> SerpentineValidator<T> {
    /// Create new validator
    pub fn new(geometry: SerpentineGeometry<T>) -> Self {
        Self { geometry }
    }

    /// Validate mixing efficiency against expected values
    ///
    /// # Validation Criteria
    ///
    /// For well-designed serpentines:
    /// - Mixing length < total channel length
    /// - Mixing fraction at outlet > 90%
    /// - Pressure drop < 100 kPa (for ~mm-scale channels)
    pub fn validate_mixing(
        &self,
        solution: &SerpentineMixingSolution<T>,
    ) -> Result<SerpentineValidationResult<T>, String> {
        let total_length = self.geometry.total_length();
        let l_mix = solution.l_mix_90;

        let achievable = l_mix < total_length * from_f64::<T>(2.0);
        let well_mixed = solution.is_well_mixed();
        let dp_reasonable = solution.pressure_drop < from_f64::<T>(1e5);
        let validation_passed = achievable && well_mixed && dp_reasonable;

        let mut error_message = None;
        if !achievable {
            error_message = Some("Mixing length exceeds 2× total channel length".to_string());
        } else if !well_mixed {
            error_message = Some(format!(
                "Mixing fraction {:.2}% < 90% at outlet",
                eunomia::NumericElement::to_f64(
                    solution.mixing_fraction_outlet * from_f64::<T>(100.0),
                )
            ));
        } else if !dp_reasonable {
            error_message = Some(format!(
                "Pressure drop {:.0} Pa > 100 kPa",
                eunomia::NumericElement::to_f64(solution.pressure_drop)
            ));
        }

        Ok(SerpentineValidationResult {
            validation_passed,
            mixing_fraction_outlet: solution.mixing_fraction_outlet,
            error_message,
        })
    }
}

/// Validation result for serpentine
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SerpentineValidationResult<T: Cfd2dScalar + Copy> {
    /// Validation passed
    pub validation_passed: bool,
    /// Measured mixing fraction at outlet
    pub mixing_fraction_outlet: T,
    /// Error message if validation failed
    pub error_message: Option<String>,
}
