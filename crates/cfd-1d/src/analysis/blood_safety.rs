//! Blood-safety analysis utilities for shear-limit screening.
//!
//! This module provides configurable limit checks to flag channel segments that
//! may exceed blood-handling shear thresholds used during device risk screening.

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Configurable blood shear limits used to flag potentially unsafe conditions.
#[derive(Debug, Clone)]
pub struct BloodShearLimits<T: RealField + Copy> {
    /// Maximum allowable wall shear stress [Pa].
    pub max_wall_shear_stress_pa: T,
    /// Optional maximum allowable wall shear rate [1/s].
    pub max_wall_shear_rate_per_s: Option<T>,
}

impl<T: RealField + Copy + FromPrimitive> BloodShearLimits<T> {
    /// Conservative default profile for FDA-oriented whole-blood design screening.
    ///
    /// FDA guidance is risk-based and does not prescribe a single universal shear
    /// cutoff for all devices; this default is intentionally conservative and should
    /// be tuned to the specific indication, exposure duration, and test evidence.
    #[must_use]
    pub fn fda_conservative_whole_blood() -> Self {
        Self {
            max_wall_shear_stress_pa: T::from_f64(150.0).unwrap_or_else(T::one),
            max_wall_shear_rate_per_s: None,
        }
    }
}

/// Shear-limit violation details for a single component.
#[derive(Debug, Clone)]
pub struct ShearLimitViolation<T: RealField + Copy> {
    /// Component identifier.
    pub component_id: String,
    /// Computed wall shear stress [Pa].
    pub wall_shear_stress_pa: T,
    /// Stress limit that was exceeded [Pa].
    pub stress_limit_pa: T,
    /// Stress exceedance ratio = wall_shear_stress_pa / stress_limit_pa.
    pub stress_exceedance_ratio: T,
    /// Computed wall shear rate [1/s], if available.
    pub wall_shear_rate_per_s: Option<T>,
    /// Shear-rate limit [1/s], if configured.
    pub shear_rate_limit_per_s: Option<T>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservative_fda_profile_has_positive_stress_limit() {
        let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood();
        assert!(limits.max_wall_shear_stress_pa > 0.0);
    }
}
