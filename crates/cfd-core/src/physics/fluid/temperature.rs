//! Temperature-dependent fluid models

use super::traits::{Fluid as FluidTrait, FluidState};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Temperature-dependent viscosity model using polynomial fit
///
/// μ(T) = `Σ(a_i` * T^i) for i = 0 to n
///
/// NOTE: Currently lacks a framework for fitting these polynomials from real fluid data.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PolynomialViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Base density [kg/m³]
    pub density_ref: T,
    /// Temperature coefficient for density [1/K]
    pub thermal_expansion: T,
    /// Polynomial coefficients for viscosity [Pa·s/K^i]
    pub viscosity_coeffs: Vec<T>,
    /// Reference temperature [K]
    pub t_ref: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + FromPrimitive + Copy> PolynomialViscosity<T> {
    /// Calculate viscosity at given temperature
    pub fn calculate_viscosity(&self, temperature: T) -> T {
        let mut viscosity = T::zero();
        let mut t_power = T::one();

        for coeff in &self.viscosity_coeffs {
            viscosity += *coeff * t_power;
            t_power *= temperature;
        }

        viscosity
    }

    /// Calculate density with thermal expansion
    pub fn calculate_density(&self, temperature: T) -> T {
        let dt = temperature - self.t_ref;
        self.density_ref * (T::one() - self.thermal_expansion * dt)
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for PolynomialViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let density = self.calculate_density(temperature);
        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }

    fn reference_temperature(&self) -> Option<T> {
        Some(self.t_ref)
    }
}

// TODO: Implement Arrhenius viscosity model (μ(T) = A * exp(B/T))
// Common for liquids where viscosity decreases with temperature

/// Andrade viscosity model
///
/// μ(T) = A * exp(B/(T - C))
/// Extended version of Arrhenius with additional parameter
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AndradeViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Pre-exponential factor A [Pa·s]
    pub a_factor: T,
    /// Temperature coefficient B [K]
    pub b_factor: T,
    /// Temperature offset C [K]
    pub c_factor: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + FromPrimitive + Copy> AndradeViscosity<T> {
    /// Calculate viscosity using Andrade model
    pub fn calculate_viscosity(&self, temperature: T) -> Result<T, Error> {
        let denominator = temperature - self.c_factor;
        if denominator <= T::zero() {
            return Err(Error::InvalidInput(format!(
                "Temperature must be greater than {}",
                self.c_factor
            )));
        }

        Ok(self.a_factor * (self.b_factor / denominator).exp())
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for AndradeViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let viscosity = self.calculate_viscosity(temperature)?;

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }
}

/// Sutherland viscosity model for gases
///
/// μ(T) = μ_ref * (T/T_ref)^(3/2) * (T_ref + S) / (T + S)
/// Common for gases (e.g., Air)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SutherlandViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³] (assumed constant for this model)
    pub density: T,
    /// Reference viscosity [Pa·s]
    pub mu_ref: T,
    /// Reference temperature [K]
    pub t_ref: T,
    /// Sutherland constant [K]
    pub s_constant: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + FromPrimitive + Copy> SutherlandViscosity<T> {
    /// Calculate viscosity using Sutherland's law
    pub fn calculate_viscosity(&self, temperature: T) -> T {
        let t_ratio = temperature / self.t_ref;
        let numerator = self.t_ref + self.s_constant;
        let denominator = temperature + self.s_constant;

        self.mu_ref * t_ratio.powf(T::from_f64(1.5).unwrap_or_else(T::one)) * numerator
            / denominator
    }
}

impl<T: RealField + FromPrimitive + Copy> FluidTrait<T> for SutherlandViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
            speed_of_sound: self.speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sutherland_viscosity_air() {
        // Sutherland parameters for Air
        // mu_ref = 1.716e-5 Pa·s at T_ref = 273.15 K, S = 110.4 K
        let air = SutherlandViscosity::<f64> {
            name: "Air".to_string(),
            density: 1.225,              // kg/m^3
            mu_ref: 1.716e-5,            // Pa·s
            t_ref: 273.15,               // K
            s_constant: 110.4,           // K
            specific_heat: 1005.0,       // J/(kg·K)
            thermal_conductivity: 0.024, // W/(m·K)
            speed_of_sound: 340.0,       // m/s
        };

        // Test at reference temperature
        let mu_ref_calc = air.calculate_viscosity(273.15);
        assert!((mu_ref_calc - 1.716e-5).abs() < 1e-10);

        // Test at higher temperature (e.g. 300 K)
        // Hand calculation:
        // T = 300
        // (300/273.15)^1.5 * (273.15 + 110.4) / (300 + 110.4) * 1.716e-5
        // 1.15016 * 383.55 / 410.4 * 1.716e-5
        // 1.15016 * 0.93457 * 1.716e-5 = 1.0749 * 1.716e-5 = 1.844e-5 roughly
        let mu_300 = air.calculate_viscosity(300.0);
        // Expected value from online calculators/tables for air at 300K is approx 1.846e-5
        assert!((mu_300 - 1.846e-5).abs() < 1e-7);
    }

    #[test]
    fn test_sutherland_negative_temperature() {
        let air = SutherlandViscosity::<f64> {
            name: "Air".to_string(),
            density: 1.225,
            mu_ref: 1.716e-5,
            t_ref: 273.15,
            s_constant: 110.4,
            specific_heat: 1005.0,
            thermal_conductivity: 0.024,
            speed_of_sound: 340.0,
        };

        assert!(air.properties_at(-10.0, 101325.0).is_err());
        assert!(air.properties_at(0.0, 101325.0).is_err());
    }

}
