//! Cavitation number calculations.

use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Cavitation number (dimensionless)
/// σ = (p - `p_v`) / (0.5 * ρ * v²)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CavitationNumber<T: FloatElement + Copy> {
    /// Reference pressure (Pa)
    pub reference_pressure: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Reference density (kg/m³)
    pub density: T,
    /// Reference velocity (m/s)
    pub velocity: T,
}

impl<T: FloatElement + Copy> CavitationNumber<T> {
    /// Calculate cavitation number
    pub fn calculate(&self) -> T {
        let half = <T as FloatElement>::from_f64(0.5);
        let dynamic_pressure = half * self.density * self.velocity * self.velocity;

        if dynamic_pressure > <T as FloatElement>::from_f64(1e-10) {
            (self.reference_pressure - self.vapor_pressure) / dynamic_pressure
        } else {
            // Large value for zero velocity.
            <T as FloatElement>::from_f64(1e10)
        }
    }

    /// Check if cavitation is likely to occur
    pub fn is_cavitating(&self) -> bool {
        let sigma = self.calculate();
        let threshold =
            <T as FloatElement>::from_f64(super::constants::CAVITATION_INCEPTION_THRESHOLD);
        sigma < threshold
    }

    /// Calculate incipient cavitation index (Thoma number)
    pub fn thoma_number(&self, head: T) -> T {
        let g = <T as FloatElement>::from_f64(9.81);
        (self.reference_pressure - self.vapor_pressure) / (self.density * g * head)
    }

    /// Calculate pressure recovery coefficient
    pub fn pressure_recovery(&self, downstream_pressure: T) -> T {
        let half = <T as FloatElement>::from_f64(0.5);
        let dynamic_pressure = half * self.density * self.velocity * self.velocity;

        if dynamic_pressure > <T as FloatElement>::from_f64(1e-10) {
            (downstream_pressure - self.reference_pressure) / dynamic_pressure
        } else {
            <T as NumericElement>::ZERO
        }
    }
}

#[cfg(test)]
mod tests {
    use super::CavitationNumber;

    fn representative_number() -> CavitationNumber<f64> {
        CavitationNumber {
            reference_pressure: 101_325.0,
            vapor_pressure: 2_339.0,
            density: 1_000.0,
            velocity: 10.0,
        }
    }

    #[test]
    fn cavitation_number_matches_definition() {
        let sigma = representative_number().calculate();
        let expected = (101_325.0 - 2_339.0) / (0.5 * 1_000.0 * 10.0 * 10.0);

        assert!((sigma - expected).abs() <= 1.0e-12);
    }

    #[test]
    fn zero_velocity_returns_large_non_cavitating_index() {
        let sigma = CavitationNumber {
            velocity: 0.0,
            ..representative_number()
        }
        .calculate();

        assert_eq!(sigma, 1.0e10);
    }

    #[test]
    fn pressure_recovery_matches_dynamic_pressure_scaling() {
        let recovery = representative_number().pressure_recovery(120_000.0);
        let expected = (120_000.0 - 101_325.0) / (0.5 * 1_000.0 * 10.0 * 10.0);

        assert!((recovery - expected).abs() <= 1.0e-12);
    }
}
