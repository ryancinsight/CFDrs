//! Flow regime classification
//!
//! Provides flow regime identification based on dimensionless numbers
//! and flow characteristics.

use eunomia::{NumericElement, RealField};
use serde::{Deserialize, Serialize};

/// Flow regime classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlowRegime {
    /// Stokes flow (Re << 1)
    Stokes,
    /// Laminar flow
    Laminar,
    /// Transitional flow
    Transitional,
    /// Turbulent flow
    Turbulent,
    /// Hypersonic flow
    Hypersonic,
}

/// Flow classifier based on dimensionless numbers
pub struct FlowClassifier;

impl FlowClassifier {
    /// Classify flow regime based on Reynolds number
    pub fn classify_by_reynolds<T: RealField>(reynolds: T) -> FlowRegime {
        let re = <T as NumericElement>::to_f64(reynolds);

        if re < 0.1 {
            FlowRegime::Stokes
        } else if re < crate::physics::constants::physics::dimensionless::reynolds::PIPE_LAMINAR_MAX
        {
            FlowRegime::Laminar
        } else if re
            < crate::physics::constants::physics::dimensionless::reynolds::PIPE_TURBULENT_MIN
        {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }

    /// Classify flow regime based on Mach number
    pub fn classify_by_mach<T: RealField>(mach: T) -> FlowRegime {
        let ma = <T as NumericElement>::to_f64(mach);

        if ma < crate::physics::constants::physics::dimensionless::mach::INCOMPRESSIBLE_LIMIT {
            // Incompressible flow - further classify by Reynolds
            FlowRegime::Laminar // Would need Reynolds for proper classification
        } else if ma < crate::physics::constants::physics::dimensionless::mach::TRANSONIC_LOWER {
            // Subsonic compressible
            FlowRegime::Turbulent // Typically turbulent at high speeds
        } else if ma < crate::physics::constants::physics::dimensionless::mach::HYPERSONIC {
            // Supersonic
            FlowRegime::Turbulent
        } else {
            // Hypersonic
            FlowRegime::Hypersonic
        }
    }

    /// Classify based on multiple dimensionless numbers
    pub fn classify<T: RealField>(reynolds: T, mach: Option<T>, _froude: Option<T>) -> FlowRegime {
        // Priority: Mach number for compressibility, then Reynolds for turbulence
        if let Some(ma) = mach {
            if <T as NumericElement>::to_f64(ma)
                >= crate::physics::constants::physics::dimensionless::mach::HYPERSONIC
            {
                return FlowRegime::Hypersonic;
            }
        }

        // Use Reynolds number for flow regime classification
        Self::classify_by_reynolds(reynolds)
    }
}

impl FlowRegime {
    /// Check if flow is viscous-dominated
    #[must_use]
    pub fn is_viscous_dominated(&self) -> bool {
        matches!(self, FlowRegime::Stokes | FlowRegime::Laminar)
    }

    /// Check if flow is inertia-dominated
    #[must_use]
    pub fn is_inertia_dominated(&self) -> bool {
        matches!(self, FlowRegime::Turbulent | FlowRegime::Hypersonic)
    }

    /// Check if flow requires turbulence modeling
    #[must_use]
    pub fn requires_turbulence_model(&self) -> bool {
        matches!(self, FlowRegime::Transitional | FlowRegime::Turbulent)
    }

    /// Get typical CFL number for this regime
    #[must_use]
    pub fn typical_cfl(&self) -> f64 {
        match self {
            FlowRegime::Stokes => 10.0, // Can use large time steps
            FlowRegime::Laminar => 1.0,
            FlowRegime::Transitional => 0.5,
            FlowRegime::Turbulent => 0.3,
            FlowRegime::Hypersonic => 0.1, // Need small time steps
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{FlowClassifier, FlowRegime};

    #[test]
    fn reynolds_classification_respects_pipe_thresholds() {
        assert_eq!(
            FlowClassifier::classify_by_reynolds(0.05_f64),
            FlowRegime::Stokes
        );
        assert_eq!(
            FlowClassifier::classify_by_reynolds(1_000.0_f64),
            FlowRegime::Laminar
        );
        assert_eq!(
            FlowClassifier::classify_by_reynolds(3_000.0_f64),
            FlowRegime::Transitional
        );
        assert_eq!(
            FlowClassifier::classify_by_reynolds(5_000.0_f64),
            FlowRegime::Turbulent
        );
    }

    #[test]
    fn mach_classification_uses_compressibility_thresholds() {
        assert_eq!(
            FlowClassifier::classify_by_mach(0.1_f32),
            FlowRegime::Laminar
        );
        assert_eq!(
            FlowClassifier::classify_by_mach(0.5_f32),
            FlowRegime::Turbulent
        );
        assert_eq!(
            FlowClassifier::classify_by_mach(1.5_f32),
            FlowRegime::Turbulent
        );
        assert_eq!(
            FlowClassifier::classify_by_mach(6.0_f32),
            FlowRegime::Hypersonic
        );
    }

    #[test]
    fn combined_classification_prioritizes_hypersonic_mach() {
        assert_eq!(
            FlowClassifier::classify(100.0_f64, Some(6.0_f64), None),
            FlowRegime::Hypersonic
        );
        assert_eq!(
            FlowClassifier::classify(100.0_f64, Some(0.2_f64), None),
            FlowRegime::Laminar
        );
    }
}
