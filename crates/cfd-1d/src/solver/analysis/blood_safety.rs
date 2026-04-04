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
    /// Optional maximum allowable Giersiepen hemolysis index [-].
    pub max_giersiepen_hi: Option<T>,
    /// Optional maximum allowable Taskin hemolysis index [-].
    pub max_taskin_hi: Option<T>,
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
            max_wall_shear_stress_pa: T::from_f64(150.0)
                .expect("Mathematical constant conversion compromised"),
            max_wall_shear_rate_per_s: None,
            max_giersiepen_hi: None,
            max_taskin_hi: None,
        }
    }

    /// Attach optional exposure-time-aware hemolysis limits to the screening profile.
    #[must_use]
    pub fn with_hemolysis_limits(
        mut self,
        max_giersiepen_hi: Option<T>,
        max_taskin_hi: Option<T>,
    ) -> Self {
        self.max_giersiepen_hi = max_giersiepen_hi;
        self.max_taskin_hi = max_taskin_hi;
        self
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

/// Exposure-time-aware hemolysis-limit violation details for a single component.
#[derive(Debug, Clone)]
pub struct HemolysisLimitViolation<T: RealField + Copy> {
    /// Component identifier.
    pub component_id: String,
    /// Computed wall shear stress [Pa].
    pub wall_shear_stress_pa: T,
    /// Exposure duration used to evaluate damage accumulation [s].
    pub exposure_time_s: T,
    /// Computed Giersiepen hemolysis index [-], if configured.
    pub giersiepen_hi: Option<T>,
    /// Configured Giersiepen limit [-], if configured.
    pub giersiepen_limit: Option<T>,
    /// Giersiepen exceedance ratio = giersiepen_hi / giersiepen_limit.
    pub giersiepen_exceedance_ratio: Option<T>,
    /// Computed Taskin hemolysis index [-], if configured.
    pub taskin_hi: Option<T>,
    /// Configured Taskin limit [-], if configured.
    pub taskin_limit: Option<T>,
    /// Taskin exceedance ratio = taskin_hi / taskin_limit.
    pub taskin_exceedance_ratio: Option<T>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservative_fda_profile_has_positive_stress_limit() {
        let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood();
        assert!(limits.max_wall_shear_stress_pa > 0.0);
        assert!(limits.max_giersiepen_hi.is_none());
        assert!(limits.max_taskin_hi.is_none());
    }

    #[test]
    fn with_hemolysis_limits_overrides_defaults() {
        let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood()
            .with_hemolysis_limits(Some(1e-4), Some(5e-4));

        assert_eq!(limits.max_giersiepen_hi, Some(1e-4));
        assert_eq!(limits.max_taskin_hi, Some(5e-4));
    }
}
