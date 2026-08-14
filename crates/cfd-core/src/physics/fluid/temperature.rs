//! Temperature-dependent fluid models

#![cfg_attr(
    test,
    expect(
        clippy::unwrap_used,
        reason = "ratchet CFDRS-UNWRAP-1: pre-existing debt"
    )
)]

use super::thermophysical::linear_density_at;
use super::traits::{Fluid as FluidTrait, FluidState};
use crate::error::Error;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, ReciprocalTemperature, SpecificHeatCapacity,
    TemperatureDifference, ThermalConductivity, ThermodynamicTemperature, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
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
    pub density_ref: MassDensity<T>,
    /// Temperature coefficient for density [1/K]
    pub thermal_expansion: ReciprocalTemperature<T>,
    /// Polynomial coefficients for viscosity [Pa·s/K^i] at the formula boundary.
    ///
    /// Each coefficient has a different temperature exponent and therefore a
    /// different derived unit; the vector remains scalar formula data rather
    /// than pretending all coefficients share one quantity dimension.
    pub viscosity_coeffs: Vec<T>,
    /// Reference temperature \[K]
    pub t_ref: ThermodynamicTemperature<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + NumericElement + Copy> PolynomialViscosity<T> {
    /// Calculate viscosity at given temperature
    pub fn calculate_viscosity(&self, temperature: ThermodynamicTemperature<T>) -> T {
        let temperature = temperature.into_base();
        let mut viscosity = <T as NumericElement>::ZERO;
        let mut t_power = <T as NumericElement>::ONE;

        for coeff in &self.viscosity_coeffs {
            viscosity += *coeff * t_power;
            t_power *= temperature;
        }

        viscosity
    }

    /// Calculate density with the shared Proteus linear temperature response.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidInput`] when the reference thermophysical state,
    /// response coefficient, evaluation temperature, or resulting density lies
    /// outside the Proteus material-property contract.
    pub fn calculate_density(&self, temperature: ThermodynamicTemperature<T>) -> Result<T, Error> {
        linear_density_at(
            self.density_ref.into_base(),
            self.thermal_expansion.into_base(),
            self.t_ref.into_base(),
            self.specific_heat.into_base(),
            self.thermal_conductivity.into_base(),
            temperature.into_base(),
        )
    }
}

impl<T: RealField + NumericElement + Copy> FluidTrait<T> for PolynomialViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let temperature = ThermodynamicTemperature::from_base(temperature);
        if temperature.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let density = self.calculate_density(temperature)?;
        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density: MassDensity::from_base(density),
            dynamic_viscosity: DynamicViscosity::from_base(viscosity),
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

    fn reference_temperature(
        &self,
    ) -> Option<aequitas::systems::si::quantities::ThermodynamicTemperature<T>> {
        Some(self.t_ref)
    }
}

/// Arrhenius viscosity model
///
/// μ(T) = A * exp(B/T)
///
/// Common for liquids where viscosity decreases with temperature.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ArrheniusViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: MassDensity<T>,
    /// Pre-exponential factor A [Pa·s]
    pub a_factor: DynamicViscosity<T>,
    /// Temperature coefficient B \[K]
    pub b_factor: TemperatureDifference<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + FloatElement + Copy> ArrheniusViscosity<T> {
    /// Calculate viscosity using Arrhenius model
    pub fn calculate_viscosity(
        &self,
        temperature: ThermodynamicTemperature<T>,
    ) -> Result<T, Error> {
        let temperature = temperature.into_base();
        if temperature <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        Ok(self.a_factor.into_base() * FloatElement::exp(self.b_factor.into_base() / temperature))
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for ArrheniusViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let temperature = ThermodynamicTemperature::from_base(temperature);
        let viscosity = self.calculate_viscosity(temperature)?;

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: DynamicViscosity::from_base(viscosity),
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

/// Andrade viscosity model
///
/// μ(T) = A * exp(B/(T - C))
/// Extended version of Arrhenius with additional parameter
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AndradeViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: MassDensity<T>,
    /// Pre-exponential factor A [Pa·s]
    pub a_factor: DynamicViscosity<T>,
    /// Temperature coefficient B \[K]
    pub b_factor: TemperatureDifference<T>,
    /// Temperature offset C \[K]
    pub c_factor: TemperatureDifference<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + FloatElement + Copy> AndradeViscosity<T> {
    /// Calculate viscosity using Andrade model
    pub fn calculate_viscosity(
        &self,
        temperature: ThermodynamicTemperature<T>,
    ) -> Result<T, Error> {
        let temperature = temperature.into_base();
        let denominator = temperature - self.c_factor.into_base();
        if denominator <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(format!(
                "Temperature must be greater than {:?}",
                self.c_factor.into_base()
            )));
        }

        Ok(self.a_factor.into_base() * FloatElement::exp(self.b_factor.into_base() / denominator))
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for AndradeViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let temperature = ThermodynamicTemperature::from_base(temperature);
        let viscosity = self.calculate_viscosity(temperature)?;

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: DynamicViscosity::from_base(viscosity),
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
/// μ(T) = `μ_ref` * (`T/T_ref)^(3/2`) * (`T_ref` + S) / (T + S)
/// Common for gases (e.g., Air)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SutherlandViscosity<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³] (assumed constant for this model)
    pub density: MassDensity<T>,
    /// Reference viscosity [Pa·s]
    pub mu_ref: DynamicViscosity<T>,
    /// Reference temperature \[K]
    pub t_ref: ThermodynamicTemperature<T>,
    /// Sutherland constant \[K]
    pub s_constant: TemperatureDifference<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
    /// Speed of sound \[m/s]
    pub speed_of_sound: Velocity<T>,
}

impl<T: RealField + FloatElement + Copy> SutherlandViscosity<T> {
    /// Calculate viscosity using Sutherland's law
    pub fn calculate_viscosity(&self, temperature: ThermodynamicTemperature<T>) -> T {
        let temperature = temperature.into_base();
        let t_ratio = temperature / self.t_ref.into_base();
        let numerator = self.t_ref.into_base() + self.s_constant.into_base();
        let denominator = temperature + self.s_constant.into_base();

        self.mu_ref.into_base()
            * FloatElement::powf(t_ratio, <T as FloatElement>::from_f64(1.5))
            * numerator
            / denominator
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for SutherlandViscosity<T> {
    fn properties_at(&self, temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        let temperature = ThermodynamicTemperature::from_base(temperature);
        if temperature.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        let viscosity = self.calculate_viscosity(temperature);

        Ok(FluidState {
            density: self.density,
            dynamic_viscosity: DynamicViscosity::from_base(viscosity),
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
            density: MassDensity::from_base(1.225),
            mu_ref: DynamicViscosity::from_base(1.716e-5),
            t_ref: ThermodynamicTemperature::from_base(273.15),
            s_constant: TemperatureDifference::from_base(110.4),
            specific_heat: SpecificHeatCapacity::from_base(1005.0),
            thermal_conductivity: ThermalConductivity::from_base(0.024),
            speed_of_sound: Velocity::from_base(340.0),
        };

        // Test at reference temperature
        let mu_ref_calc = air.calculate_viscosity(ThermodynamicTemperature::from_base(273.15));
        assert!((mu_ref_calc - 1.716e-5).abs() < 1e-10);

        // Test at higher temperature (e.g. 300 K)
        // Hand calculation:
        // T = 300
        // (300/273.15)^1.5 * (273.15 + 110.4) / (300 + 110.4) * 1.716e-5
        // 1.15016 * 383.55 / 410.4 * 1.716e-5
        // 1.15016 * 0.93457 * 1.716e-5 = 1.0749 * 1.716e-5 = 1.844e-5 roughly
        let mu_300 = air.calculate_viscosity(ThermodynamicTemperature::from_base(300.0));
        let expected_mu_300 =
            1.716e-5 * f64::powf(300.0 / 273.15, 1.5) * (273.15 + 110.4) / (300.0 + 110.4);
        assert!((mu_300 - expected_mu_300).abs() <= 16.0 * f64::EPSILON * expected_mu_300);
    }

    #[test]
    fn test_sutherland_negative_temperature() {
        let air = SutherlandViscosity::<f64> {
            name: "Air".to_string(),
            density: MassDensity::from_base(1.225),
            mu_ref: DynamicViscosity::from_base(1.716e-5),
            t_ref: ThermodynamicTemperature::from_base(273.15),
            s_constant: TemperatureDifference::from_base(110.4),
            specific_heat: SpecificHeatCapacity::from_base(1005.0),
            thermal_conductivity: ThermalConductivity::from_base(0.024),
            speed_of_sound: Velocity::from_base(340.0),
        };

        assert!(air.properties_at(-10.0, 101325.0).is_err());
        assert!(air.properties_at(0.0, 101325.0).is_err());
    }

    #[test]
    fn test_arrhenius_viscosity() {
        let fluid = ArrheniusViscosity::<f64> {
            name: "Test".to_string(),
            density: MassDensity::from_base(1000.0),
            a_factor: DynamicViscosity::from_base(1.0),
            b_factor: TemperatureDifference::from_base(2.0),
            specific_heat: SpecificHeatCapacity::from_base(1000.0),
            thermal_conductivity: ThermalConductivity::from_base(1.0),
            speed_of_sound: Velocity::from_base(1500.0),
        };

        let mu = fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(2.0))
            .unwrap();
        assert!((mu - std::f64::consts::E).abs() < 1e-12);
        assert!(fluid.properties_at(2.0, 1.0).is_ok());
    }

    #[test]
    fn test_arrhenius_invalid_temperature() {
        let fluid = ArrheniusViscosity::<f64> {
            name: "Test".to_string(),
            density: MassDensity::from_base(1000.0),
            a_factor: DynamicViscosity::from_base(1.0),
            b_factor: TemperatureDifference::from_base(2.0),
            specific_heat: SpecificHeatCapacity::from_base(1000.0),
            thermal_conductivity: ThermalConductivity::from_base(1.0),
            speed_of_sound: Velocity::from_base(1500.0),
        };

        assert!(fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(0.0))
            .is_err());
        assert!(fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(-1.0))
            .is_err());
    }

    #[test]
    fn test_polynomial_viscosity_and_density_are_value_semantic() {
        let fluid = PolynomialViscosity::<f64> {
            name: "Polynomial".to_string(),
            density_ref: MassDensity::from_base(1000.0),
            thermal_expansion: ReciprocalTemperature::from_base(2.0e-4),
            viscosity_coeffs: vec![1.0, -0.01, 0.0001],
            t_ref: ThermodynamicTemperature::from_base(300.0),
            specific_heat: SpecificHeatCapacity::from_base(4180.0),
            thermal_conductivity: ThermalConductivity::from_base(0.6),
            speed_of_sound: Velocity::from_base(1480.0),
        };

        let viscosity = fluid.calculate_viscosity(ThermodynamicTemperature::from_base(310.0));
        let expected_viscosity = 1.0 - 0.01 * 310.0 + 0.0001 * 310.0 * 310.0;
        assert!((viscosity - expected_viscosity).abs() <= 16.0 * f64::EPSILON);

        let density = fluid
            .calculate_density(ThermodynamicTemperature::from_base(310.0))
            .expect("finite positive thermophysical state");
        let expected_density = 1000.0 * (1.0 - 2.0e-4 * 10.0);
        // Proteus evaluates the response with an FMA. Four native roundings
        // bound that path, this independent expression, and density scaling.
        let rounding = 4.0 * f64::EPSILON * expected_density.abs();
        assert!((density - expected_density).abs() <= rounding);
    }

    #[test]
    fn polynomial_density_rejects_a_negative_proteus_state() {
        let fluid = PolynomialViscosity::<f64> {
            name: "Invalid response".to_string(),
            density_ref: MassDensity::from_base(1000.0),
            thermal_expansion: ReciprocalTemperature::from_base(0.2),
            viscosity_coeffs: vec![1.0],
            t_ref: ThermodynamicTemperature::from_base(300.0),
            specific_heat: SpecificHeatCapacity::from_base(4180.0),
            thermal_conductivity: ThermalConductivity::from_base(0.6),
            speed_of_sound: Velocity::from_base(1480.0),
        };

        let error = fluid
            .calculate_density(ThermodynamicTemperature::from_base(310.0))
            .expect_err("the linear response produces negative density");
        match error {
            Error::InvalidInput(message) => {
                assert!(message.contains("MassDensity"));
                assert!(message.contains("FiniteNonNegative"));
            }
            other => panic!("unexpected error variant: {other}"),
        }
    }

    #[test]
    fn test_andrade_viscosity_and_domain_are_value_semantic() {
        let fluid = AndradeViscosity::<f64> {
            name: "Andrade".to_string(),
            density: MassDensity::from_base(1000.0),
            a_factor: DynamicViscosity::from_base(0.25),
            b_factor: TemperatureDifference::from_base(10.0),
            c_factor: TemperatureDifference::from_base(2.0),
            specific_heat: SpecificHeatCapacity::from_base(1000.0),
            thermal_conductivity: ThermalConductivity::from_base(1.0),
            speed_of_sound: Velocity::from_base(1500.0),
        };

        let viscosity = fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(7.0))
            .unwrap();
        let expected_viscosity = 0.25 * f64::exp(10.0 / (7.0 - 2.0));
        assert!((viscosity - expected_viscosity).abs() <= 8.0 * f64::EPSILON * expected_viscosity);

        assert!(fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(2.0))
            .is_err());
        assert!(fluid
            .calculate_viscosity(ThermodynamicTemperature::from_base(1.0))
            .is_err());
    }
}
