//! Blood-safety analysis utilities for shear-limit screening.
//!
//! This module provides configurable limit checks to flag channel segments that
//! may exceed blood-handling shear thresholds used during device risk screening.

use crate::scalar::Cfd1dScalar;
use aequitas::systems::si::quantities::{Dimensionless, Pressure, ReciprocalTime, Time};
use cfd_core::conversion::SafeFromF64;

/// Configurable blood shear limits used to flag potentially unsafe conditions.
#[derive(Debug, Clone)]
pub struct BloodShearLimits<T: Cfd1dScalar + Copy> {
    /// Maximum allowable wall shear stress \[Pa].
    pub max_wall_shear_stress_pa: Pressure<T>,
    /// Optional maximum allowable wall shear rate [1/s].
    pub max_wall_shear_rate_per_s: Option<ReciprocalTime<T>>,
    /// Optional maximum allowable Giersiepen hemolysis index [-].
    pub max_giersiepen_hi: Option<Dimensionless<T>>,
    /// Optional maximum allowable Taskin hemolysis index [-].
    pub max_taskin_hi: Option<Dimensionless<T>>,
}

impl<T: Cfd1dScalar + Copy + SafeFromF64> BloodShearLimits<T> {
    /// Conservative default profile for FDA-oriented whole-blood design screening.
    ///
    /// FDA guidance is risk-based and does not prescribe a single universal shear
    /// cutoff for all devices; this default is intentionally conservative and should
    /// be tuned to the specific indication, exposure duration, and test evidence.
    #[must_use]
    pub fn fda_conservative_whole_blood() -> Self {
        Self {
            max_wall_shear_stress_pa: Pressure::from_base(T::from_f64_or_one(150.0)),
            max_wall_shear_rate_per_s: None,
            max_giersiepen_hi: None,
            max_taskin_hi: None,
        }
    }

    /// Attach optional exposure-time-aware hemolysis limits to the screening profile.
    #[must_use]
    pub fn with_hemolysis_limits(
        mut self,
        max_giersiepen_hi: Option<Dimensionless<T>>,
        max_taskin_hi: Option<Dimensionless<T>>,
    ) -> Self {
        self.max_giersiepen_hi = max_giersiepen_hi;
        self.max_taskin_hi = max_taskin_hi;
        self
    }
}

/// Shear-limit violation details for a single component.
#[derive(Debug, Clone)]
pub struct ShearLimitViolation<T: Cfd1dScalar + Copy> {
    /// Component identifier.
    pub component_id: String,
    /// Computed wall shear stress \[Pa].
    pub wall_shear_stress_pa: Pressure<T>,
    /// Stress limit that was exceeded \[Pa].
    pub stress_limit_pa: Pressure<T>,
    /// Stress exceedance ratio = wall_shear_stress_pa / stress_limit_pa.
    pub stress_exceedance_ratio: Dimensionless<T>,
    /// Computed wall shear rate [1/s], if available.
    pub wall_shear_rate_per_s: Option<ReciprocalTime<T>>,
    /// Shear-rate limit [1/s], if configured.
    pub shear_rate_limit_per_s: Option<ReciprocalTime<T>>,
}

/// Exposure-time-aware hemolysis-limit violation details for a single component.
#[derive(Debug, Clone)]
pub struct HemolysisLimitViolation<T: Cfd1dScalar + Copy> {
    /// Component identifier.
    pub component_id: String,
    /// Computed wall shear stress \[Pa].
    pub wall_shear_stress_pa: Pressure<T>,
    /// Exposure duration used to evaluate damage accumulation \[s].
    pub exposure_time_s: Time<T>,
    /// Computed Giersiepen hemolysis index [-], if configured.
    pub giersiepen_hi: Option<Dimensionless<T>>,
    /// Configured Giersiepen limit [-], if configured.
    pub giersiepen_limit: Option<Dimensionless<T>>,
    /// Giersiepen exceedance ratio = giersiepen_hi / giersiepen_limit.
    pub giersiepen_exceedance_ratio: Option<Dimensionless<T>>,
    /// Computed Taskin hemolysis index [-], if configured.
    pub taskin_hi: Option<Dimensionless<T>>,
    /// Configured Taskin limit [-], if configured.
    pub taskin_limit: Option<Dimensionless<T>>,
    /// Taskin exceedance ratio = taskin_hi / taskin_limit.
    pub taskin_exceedance_ratio: Option<Dimensionless<T>>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn conservative_fda_profile_has_positive_stress_limit() {
        let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood();
        assert!(limits.max_wall_shear_stress_pa.into_base() > 0.0);
        assert!(limits.max_giersiepen_hi.is_none());
        assert!(limits.max_taskin_hi.is_none());
    }

    #[test]
    fn with_hemolysis_limits_overrides_defaults() {
        let limits = BloodShearLimits::<f64>::fda_conservative_whole_blood()
            .with_hemolysis_limits(
                Some(Dimensionless::from_base(1e-4)),
                Some(Dimensionless::from_base(5e-4)),
            );

        assert_eq!(limits.max_giersiepen_hi, Some(Dimensionless::from_base(1e-4)));
        assert_eq!(limits.max_taskin_hi, Some(Dimensionless::from_base(5e-4)));
    }
}
