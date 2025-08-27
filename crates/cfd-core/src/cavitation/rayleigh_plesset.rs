//! Rayleigh-Plesset bubble dynamics model.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Rayleigh-Plesset bubble dynamics model
///
/// Describes the growth and collapse of a spherical bubble in an infinite liquid
/// based on Rayleigh (1917) and Plesset (1949).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct RayleighPlesset<T: RealField + Copy> {
    /// Initial bubble radius (m)
    pub initial_radius: T,
    /// Liquid density (kg/m³)
    pub liquid_density: T,
    /// Liquid viscosity (Pa·s)
    pub liquid_viscosity: T,
    /// Surface tension (N/m)
    pub surface_tension: T,
    /// Vapor pressure (Pa)
    pub vapor_pressure: T,
    /// Polytropic index for gas content
    pub polytropic_index: T,
}

impl<T: RealField + Copy + FromPrimitive> RayleighPlesset<T> {
    /// Calculate bubble wall acceleration (d²R/dt²)
    ///
    /// Uses the full Rayleigh-Plesset equation including viscous and surface tension effects
    pub fn bubble_acceleration(&self, radius: T, velocity: T, ambient_pressure: T) -> T {
        // Avoid division by zero
        if radius < T::from_f64(super::constants::MIN_BUBBLE_RADIUS).unwrap_or_else(|| T::zero()) {
            return T::zero();
        }

        // Pressure inside bubble (including surface tension)
        let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());
        let four = two * two;
        let three_halves = T::from_f64(1.5).unwrap_or_else(|| T::one());

        let p_bubble = self.vapor_pressure + two * self.surface_tension / radius;

        // Pressure difference driving bubble motion
        let pressure_diff = p_bubble - ambient_pressure;

        // Rayleigh-Plesset equation terms
        let term1 = pressure_diff / (self.liquid_density * radius);
        let term2 = -three_halves * velocity * velocity / radius;
        let term3 =
            -four * self.liquid_viscosity * velocity / (self.liquid_density * radius * radius);
        let term4 = -two * self.surface_tension / (self.liquid_density * radius * radius);

        term1 + term2 + term3 + term4
    }

    /// Calculate bubble growth rate for simplified case (neglecting viscosity)
    pub fn growth_rate_inviscid(&self, _radius: T, ambient_pressure: T) -> T {
        let pressure_difference = self.vapor_pressure - ambient_pressure;

        if pressure_difference <= T::zero() {
            return T::zero(); // No growth if ambient pressure exceeds vapor pressure
        }

        // Simplified Rayleigh equation: dR/dt = sqrt(2/3 * Δp/ρ)
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap_or_else(|| T::one());
        two_thirds * (pressure_difference / self.liquid_density).sqrt()
    }

    /// Calculate collapse time from Rayleigh collapse formula
    pub fn collapse_time(&self, initial_radius: T, pressure_difference: T) -> T {
        // Rayleigh collapse time: t_c = 0.915 * R_0 * sqrt(ρ/Δp)
        let coefficient = T::from_f64(0.915).unwrap_or_else(|| T::one());

        if pressure_difference > T::zero() {
            coefficient * initial_radius * (self.liquid_density / pressure_difference).sqrt()
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Very large time for no pressure difference
        }
    }

    /// Calculate maximum bubble radius during growth (Rayleigh-Plesset)
    pub fn maximum_radius(&self, pressure_ratio: T) -> T {
        // R_max/R_0 = (p_∞/p_v)^(1/3γ) for isothermal growth
        let exponent =
            T::one() / (T::from_f64(3.0).unwrap_or_else(|| T::one()) * self.polytropic_index);
        self.initial_radius * pressure_ratio.powf(exponent)
    }

    /// Calculate bubble natural frequency
    pub fn natural_frequency(&self, radius: T, ambient_pressure: T) -> T {
        if radius <= T::zero() {
            return T::zero();
        }

        // ω_0 = 1/R * sqrt(3γ(p_0 - p_v)/ρ)
        let three = T::from_f64(3.0).unwrap_or_else(|| T::one() + T::one() + T::one());
        let pressure_term =
            three * self.polytropic_index * (ambient_pressure - self.vapor_pressure);

        if pressure_term > T::zero() {
            (pressure_term / self.liquid_density).sqrt() / radius
        } else {
            T::zero()
        }
    }

    /// Calculate critical Blake radius for unstable growth
    pub fn blake_critical_radius(&self, ambient_pressure: T) -> T {
        // R_c = 0.85 * (2σ/(p_∞ - p_v))
        let coefficient =
            T::from_f64(super::constants::BLAKE_CRITICAL_COEFFICIENT).unwrap_or_else(|| T::one());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());

        let pressure_diff = ambient_pressure - self.vapor_pressure;
        if pressure_diff > T::zero() {
            coefficient * two * self.surface_tension / pressure_diff
        } else {
            T::from_f64(1e10).unwrap_or_else(|| T::one()) // Very large radius
        }
    }
}
