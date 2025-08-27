//! Cavitation number calculations.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Cavitation number (dimensionless)
/// σ = (p - `p_v`) / (0.5 * ρ * v²)
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct CavitationNumber<T: RealField + Copy> {
    /// Reference pressure (Pa)
    pub reference_pressure: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Reference density (kg/m³)
    pub density: T,
    /// Reference velocity (m/s)
    pub velocity: T,
}

impl<T: RealField + Copy + FromPrimitive> CavitationNumber<T> {
    /// Calculate cavitation number
    pub fn calculate(&self) -> T {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * self.density * self.velocity * self.velocity;

        if dynamic_pressure > T::from_f64(1e-10).unwrap_or_else(|| T::zero()) {
            (self.reference_pressure - self.vapor_pressure) / dynamic_pressure
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Large value for zero velocity
        }
    }

    /// Check if cavitation is likely to occur
    pub fn is_cavitating(&self) -> bool {
        let sigma = self.calculate();
        let threshold = T::from_f64(super::constants::CAVITATION_INCEPTION_THRESHOLD)
            .unwrap_or_else(|| T::one());
        sigma < threshold
    }

    /// Calculate incipient cavitation index (Thoma number)
    pub fn thoma_number(&self, head: T) -> T {
        let g = T::from_f64(9.81).unwrap_or_else(|| T::from_f64(10.0).unwrap_or_else(|| T::one()));
        (self.reference_pressure - self.vapor_pressure) / (self.density * g * head)
    }

    /// Calculate pressure recovery coefficient
    pub fn pressure_recovery(&self, downstream_pressure: T) -> T {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * self.density * self.velocity * self.velocity;

        if dynamic_pressure > T::from_f64(1e-10).unwrap_or_else(|| T::zero()) {
            (downstream_pressure - self.reference_pressure) / dynamic_pressure
        } else {
            T::zero()
        }
    }
}
