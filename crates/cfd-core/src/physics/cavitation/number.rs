//! Cavitation number calculations.

use aequitas::systems::si::quantities::{Dimensionless, Length, MassDensity, Pressure, Velocity};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Cavitation number (dimensionless)
/// σ = (p - `p_v`) / (0.5 * ρ * v²)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CavitationNumber<T: FloatElement + Copy> {
    /// Reference pressure (Pa)
    pub reference_pressure: Pressure<T>,
    /// Vapor pressure (Pa)
    pub vapor_pressure: Pressure<T>,
    /// Reference density (kg/m³)
    pub density: MassDensity<T>,
    /// Reference velocity (m/s)
    pub velocity: Velocity<T>,
}

impl<T: FloatElement + Copy> CavitationNumber<T> {
    /// Calculate cavitation number
    pub fn calculate(&self) -> Dimensionless<T> {
        let half = <T as FloatElement>::from_f64(0.5);
        let density = self.density.into_base();
        let velocity = self.velocity.into_base();
        let dynamic_pressure = half * density * velocity * velocity;

        if dynamic_pressure > <T as FloatElement>::from_f64(1e-10) {
            Dimensionless::from_base(
                (self.reference_pressure.into_base() - self.vapor_pressure.into_base())
                    / dynamic_pressure,
            )
        } else {
            // Large value for zero velocity.
            Dimensionless::from_base(<T as FloatElement>::from_f64(1e10))
        }
    }

    /// Check if cavitation is likely to occur
    pub fn is_cavitating(&self) -> bool {
        let sigma = self.calculate().into_base();
        let threshold =
            <T as FloatElement>::from_f64(super::constants::CAVITATION_INCEPTION_THRESHOLD);
        sigma < threshold
    }

    /// Calculate incipient cavitation index (Thoma number)
    pub fn thoma_number(&self, head: Length<T>) -> Dimensionless<T> {
        let g = <T as FloatElement>::from_f64(9.81);
        Dimensionless::from_base(
            (self.reference_pressure.into_base() - self.vapor_pressure.into_base())
                / (self.density.into_base() * g * head.into_base()),
        )
    }

    /// Calculate pressure recovery coefficient
    pub fn pressure_recovery(&self, downstream_pressure: Pressure<T>) -> Dimensionless<T> {
        let half = <T as FloatElement>::from_f64(0.5);
        let dynamic_pressure =
            half * self.density.into_base() * self.velocity.into_base() * self.velocity.into_base();

        if dynamic_pressure > <T as FloatElement>::from_f64(1e-10) {
            Dimensionless::from_base(
                (downstream_pressure.into_base() - self.reference_pressure.into_base())
                    / dynamic_pressure,
            )
        } else {
            Dimensionless::from_base(<T as NumericElement>::ZERO)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::CavitationNumber;
    use aequitas::systems::si::quantities::{MassDensity, Pressure, Velocity};

    fn representative_number() -> CavitationNumber<f64> {
        CavitationNumber {
            reference_pressure: Pressure::from_base(101_325.0),
            vapor_pressure: Pressure::from_base(2_339.0),
            density: MassDensity::from_base(1_000.0),
            velocity: Velocity::from_base(10.0),
        }
    }

    #[test]
    fn cavitation_number_matches_definition() {
        let sigma = representative_number().calculate().into_base();
        let expected = (101_325.0 - 2_339.0) / (0.5 * 1_000.0 * 10.0 * 10.0);

        assert!((sigma - expected).abs() <= 1.0e-12);
    }

    #[test]
    fn zero_velocity_returns_large_non_cavitating_index() {
        let sigma = CavitationNumber {
            velocity: Velocity::from_base(0.0),
            ..representative_number()
        }
        .calculate()
        .into_base();

        assert_eq!(sigma, 1.0e10);
    }

    #[test]
    fn pressure_recovery_matches_dynamic_pressure_scaling() {
        let recovery = representative_number()
            .pressure_recovery(Pressure::from_base(120_000.0))
            .into_base();
        let expected = (120_000.0 - 101_325.0) / (0.5 * 1_000.0 * 10.0 * 10.0);

        assert!((recovery - expected).abs() <= 1.0e-12);
    }
}
