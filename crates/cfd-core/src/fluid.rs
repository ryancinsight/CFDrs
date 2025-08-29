//! Fluid properties and models with trait-based extensibility.

use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Trait defining the interface for fluid models
///
/// This trait allows for different fluid models to be implemented
/// while providing a consistent interface to the solvers.
pub trait FluidModel<T: RealField>: Send + Sync {
    /// Get fluid density [kg/m³] at given conditions
    fn density(&self, temperature: T, pressure: T) -> Result<T, Error>;

    /// Get dynamic viscosity [Pa·s] at given conditions
    fn dynamic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error>;

    /// Get kinematic viscosity [m²/s] at given conditions
    fn kinematic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error>
    where
        T: Copy,
    {
        let rho = self.density(temperature, pressure)?;
        let mu = self.dynamic_viscosity(temperature, pressure)?;
        if rho <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        Ok(mu / rho)
    }

    /// Get specific heat capacity [J/(kg·K)] at given conditions
    fn specific_heat(&self, temperature: T, pressure: T) -> Result<T, Error>;

    /// Get thermal conductivity [W/(m·K)] at given conditions
    fn thermal_conductivity(&self, temperature: T, pressure: T) -> Result<T, Error>;

    /// Get Prandtl number (dimensionless) at given conditions
    fn prandtl_number(&self, temperature: T, pressure: T) -> Result<T, Error>
    where
        T: Copy,
    {
        let cp = self.specific_heat(temperature, pressure)?;
        let mu = self.dynamic_viscosity(temperature, pressure)?;
        let k = self.thermal_conductivity(temperature, pressure)?;
        if k <= T::zero() {
            return Err(Error::InvalidInput(
                "Thermal conductivity must be positive".to_string(),
            ));
        }
        Ok(mu * cp / k)
    }

    /// Calculate Reynolds number for given flow conditions
    fn reynolds_number(
        &self,
        velocity: T,
        characteristic_length: T,
        temperature: T,
        pressure: T,
    ) -> Result<T, Error>
    where
        T: Copy,
    {
        let rho = self.density(temperature, pressure)?;
        let mu = self.dynamic_viscosity(temperature, pressure)?;
        if mu <= T::zero() {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        Ok(rho * velocity * characteristic_length / mu)
    }

    /// Get a descriptive name for this fluid model
    fn name(&self) -> &str;
}

/// Constant property fluid model (incompressible, Newtonian)
///
/// This model assumes fluid properties are independent of temperature and pressure.
/// Suitable for isothermal, incompressible flow simulations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConstantPropertyFluid<T: RealField + Copy> {
    /// Descriptive name of the fluid
    name: String,
    /// Constant density [kg/m³]
    density: T,
    /// Constant dynamic viscosity [Pa·s]
    viscosity: T,
    /// Constant specific heat capacity [J/(kg·K)]
    specific_heat: T,
    /// Constant thermal conductivity [W/(m·K)]
    thermal_conductivity: T,
}

impl<T: RealField + Copy> ConstantPropertyFluid<T> {
    /// Create a new constant property fluid
    pub fn new(
        name: String,
        density: T,
        viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
    ) -> Self {
        Self {
            name,
            density,
            viscosity,
            specific_heat,
            thermal_conductivity,
        }
    }

    /// Create water at 20°C and 1 atm
    pub fn water_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Water at 20°C, 1 atm".to_string(),
            density: T::from_f64(998.2)
                .expect("Failed to represent water density (998.2 kg/m³) in numeric type T"),
            viscosity: T::from_f64(1.002e-3)
                .expect("Failed to represent water viscosity (1.002e-3 Pa·s) in numeric type T"),
            specific_heat: T::from_f64(4182.0).expect(
                "Failed to represent water specific heat (4182 J/(kg·K)) in numeric type T",
            ),
            thermal_conductivity: T::from_f64(0.598).expect(
                "Failed to represent water thermal conductivity (0.598 W/(m·K)) in numeric type T",
            ),
        }
    }

    /// Create air at 20°C and 1 atm
    pub fn air_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Air at 20°C, 1 atm".to_string(),
            density: T::from_f64(1.204)
                .expect("Failed to represent air density (1.204 kg/m³) in numeric type T"),
            viscosity: T::from_f64(1.82e-5)
                .expect("Failed to represent air viscosity (1.82e-5 Pa·s) in numeric type T"),
            specific_heat: T::from_f64(1005.0)
                .expect("Failed to represent air specific heat (1005 J/(kg·K)) in numeric type T"),
            thermal_conductivity: T::from_f64(0.0257).expect(
                "Failed to represent air thermal conductivity (0.0257 W/(m·K)) in numeric type T",
            ),
        }
    }

    /// Get the constant density value
    pub fn density(&self) -> T {
        self.density
    }

    /// Get the constant dynamic viscosity value
    pub fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Get the kinematic viscosity
    pub fn kinematic_viscosity(&self) -> T {
        self.viscosity / self.density
    }
}

impl<T: RealField + Copy> FluidModel<T> for ConstantPropertyFluid<T> {
    fn density(&self, _temperature: T, _pressure: T) -> Result<T, Error> {
        Ok(self.density)
    }

    fn dynamic_viscosity(&self, _temperature: T, _pressure: T) -> Result<T, Error> {
        Ok(self.viscosity)
    }

    fn specific_heat(&self, _temperature: T, _pressure: T) -> Result<T, Error> {
        Ok(self.specific_heat)
    }

    fn thermal_conductivity(&self, _temperature: T, _pressure: T) -> Result<T, Error> {
        Ok(self.thermal_conductivity)
    }

    fn name(&self) -> &str {
        &self.name
    }
}

/// Ideal gas model
///
/// Models a gas using the ideal gas law: ρ = P·M/(R·T)
/// where M is molar mass and R is the universal gas constant.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct IdealGas<T: RealField + Copy> {
    /// Descriptive name
    name: String,
    /// Molar mass [kg/mol]
    molar_mass: T,
    /// Specific heat at constant pressure [J/(kg·K)]
    specific_heat_cp: T,
    /// Specific heat at constant volume [J/(kg·K)]
    specific_heat_cv: T,
    /// Reference viscosity [Pa·s] at reference temperature
    viscosity_ref: T,
    /// Reference temperature [K] for viscosity
    temperature_ref: T,
    /// Sutherland constant [K] for viscosity calculation
    sutherland_constant: T,
    /// Thermal conductivity coefficient
    thermal_conductivity_coeff: T,
}

impl<T: RealField + Copy + FromPrimitive> IdealGas<T> {
    /// Universal gas constant [J/(mol·K)]
    const R_UNIVERSAL: f64 = 8.314462618;

    /// Create a new ideal gas model
    pub fn new(
        name: String,
        molar_mass: T,
        specific_heat_cp: T,
        specific_heat_cv: T,
        viscosity_ref: T,
        temperature_ref: T,
        sutherland_constant: T,
        thermal_conductivity_coeff: T,
    ) -> Self {
        Self {
            name,
            molar_mass,
            specific_heat_cp,
            specific_heat_cv,
            viscosity_ref,
            temperature_ref,
            sutherland_constant,
            thermal_conductivity_coeff,
        }
    }

    /// Create an ideal gas model for dry air
    pub fn air() -> Self {
        Self {
            name: "Dry Air (Ideal Gas)".to_string(),
            molar_mass: T::from_f64(0.02897)
                .expect("Failed to represent air molar mass (0.02897 kg/mol)"),
            specific_heat_cp: T::from_f64(1005.0)
                .expect("Failed to represent air Cp (1005 J/(kg·K))"),
            specific_heat_cv: T::from_f64(718.0)
                .expect("Failed to represent air Cv (718 J/(kg·K))"),
            viscosity_ref: T::from_f64(1.716e-5)
                .expect("Failed to represent reference viscosity (1.716e-5 Pa·s)"),
            temperature_ref: T::from_f64(273.15)
                .expect("Failed to represent reference temperature (273.15 K)"),
            sutherland_constant: T::from_f64(110.4)
                .expect("Failed to represent Sutherland constant (110.4 K)"),
            thermal_conductivity_coeff: T::from_f64(2.624e-3)
                .expect("Failed to represent thermal conductivity coefficient"),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> FluidModel<T> for IdealGas<T> {
    fn density(&self, temperature: T, pressure: T) -> Result<T, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        let r_universal =
            T::from_f64(Self::R_UNIVERSAL).expect("Failed to represent universal gas constant");
        let r_specific = r_universal / self.molar_mass;
        Ok(pressure / (r_specific * temperature))
    }

    fn dynamic_viscosity(&self, temperature: T, _pressure: T) -> Result<T, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        // Sutherland's formula for viscosity
        let t_ratio = temperature / self.temperature_ref;
        let t_ratio_3_2 = t_ratio.powf(T::from_f64(1.5).expect("Failed to represent 1.5"));
        let numerator = self.temperature_ref + self.sutherland_constant;
        let denominator = temperature + self.sutherland_constant;
        Ok(self.viscosity_ref * t_ratio_3_2 * numerator / denominator)
    }

    fn specific_heat(&self, _temperature: T, _pressure: T) -> Result<T, Error> {
        // For ideal gas, Cp is constant
        Ok(self.specific_heat_cp)
    }

    fn thermal_conductivity(&self, temperature: T, _pressure: T) -> Result<T, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        // Simple temperature-dependent model
        let t_ratio = temperature / self.temperature_ref;
        let t_ratio_sqrt = t_ratio.sqrt();
        Ok(self.thermal_conductivity_coeff * t_ratio_sqrt)
    }

    fn name(&self) -> &str {
        &self.name
    }
}

/// Legacy Fluid struct for backward compatibility
///
/// DEPRECATED: Use ConstantPropertyFluid instead.
/// This struct is maintained only for backward compatibility and will be removed in v2.0.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[deprecated(
    since = "1.20.0",
    note = "Use ConstantPropertyFluid with FluidModel trait instead"
)]
pub struct Fluid<T: RealField + Copy> {
    /// Fluid name
    pub name: String,
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: Option<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: Option<T>,
}

#[allow(deprecated)]
impl<T: RealField + Copy> Fluid<T> {
    /// Create a fluid with required properties
    pub fn create(name: String, density: T, viscosity: T) -> Self {
        Self {
            name,
            density,
            viscosity,
            specific_heat: None,
            thermal_conductivity: None,
        }
    }

    /// Get kinematic viscosity
    pub fn kinematic_viscosity(&self) -> T {
        self.viscosity / self.density
    }

    /// Get dynamic viscosity
    pub fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Get Prandtl number (if thermal properties are set)
    pub fn prandtl_number(&self) -> Option<T> {
        match (self.specific_heat, self.thermal_conductivity) {
            (Some(cp), Some(k)) if k > T::zero() => Some(self.viscosity * cp / k),
            _ => None,
        }
    }

    /// Create water at 20°C
    pub fn water_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Water at 20°C".to_string(),
            density: T::from_f64(998.2)
                .expect("Failed to represent water density (998.2 kg/m³) in numeric type T"),
            viscosity: T::from_f64(1.002e-3)
                .expect("Failed to represent water viscosity (1.002e-3 Pa·s) in numeric type T"),
            specific_heat: Some(T::from_f64(4182.0).expect(
                "Failed to represent water specific heat (4182 J/(kg·K)) in numeric type T",
            )),
            thermal_conductivity: Some(T::from_f64(0.598).expect(
                "Failed to represent water thermal conductivity (0.598 W/(m·K)) in numeric type T",
            )),
        }
    }

    /// Create air at 20°C
    pub fn air_20c() -> Self
    where
        T: FromPrimitive,
    {
        Self {
            name: "Air at 20°C".to_string(),
            density: T::from_f64(1.204)
                .expect("Failed to represent air density (1.204 kg/m³) in numeric type T"),
            viscosity: T::from_f64(1.82e-5)
                .expect("Failed to represent air viscosity (1.82e-5 Pa·s) in numeric type T"),
            specific_heat: Some(
                T::from_f64(1005.0).expect(
                    "Failed to represent air specific heat (1005 J/(kg·K)) in numeric type T",
                ),
            ),
            thermal_conductivity: Some(T::from_f64(0.0257).expect(
                "Failed to represent air thermal conductivity (0.0257 W/(m·K)) in numeric type T",
            )),
        }
    }

    /// Calculate Reynolds number
    pub fn reynolds_number(&self, velocity: T, characteristic_length: T) -> T {
        (self.density * velocity * characteristic_length) / self.viscosity
    }

    /// Set specific heat
    pub fn with_specific_heat(mut self, cp: T) -> Self {
        self.specific_heat = Some(cp);
        self
    }

    /// Set thermal conductivity
    pub fn with_thermal_conductivity(mut self, k: T) -> Self {
        self.thermal_conductivity = Some(k);
        self
    }
}
