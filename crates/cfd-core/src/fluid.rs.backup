//! Fluid properties and models with trait-based extensibility.

use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// A block of computed fluid properties at a single state point.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FluidProperties<T: RealField + Copy> {
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub dynamic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
}

impl<T: RealField + Copy> FluidProperties<T> {
    /// Calculate kinematic viscosity [m²/s] from base properties.
    pub fn kinematic_viscosity(&self) -> Result<T, Error> {
        if self.density <= T::zero() {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        Ok(self.dynamic_viscosity / self.density)
    }

    /// Calculate Prandtl number from base properties.
    pub fn prandtl_number(&self) -> Result<T, Error> {
        if self.thermal_conductivity <= T::zero() {
            return Err(Error::InvalidInput(
                "Thermal conductivity must be positive".to_string(),
            ));
        }
        Ok(self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity)
    }

    /// Calculate Reynolds number for given flow conditions.
    pub fn reynolds_number(&self, velocity: T, characteristic_length: T) -> Result<T, Error> {
        if self.dynamic_viscosity <= T::zero() {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        Ok(self.density * velocity * characteristic_length / self.dynamic_viscosity)
    }
}

/// Trait defining the interface for fluid models
///
/// This trait allows for different fluid models to be implemented
/// while providing a consistent interface to the solvers.
pub trait FluidModel<T: RealField + Copy>: Send + Sync {
    /// Calculate a block of all fluid properties at the given conditions.
    fn properties(&self, temperature: T, pressure: T) -> Result<FluidProperties<T>, Error>;

    /// Get a descriptive name for this fluid model.
    fn name(&self) -> &str;

    // Convenience methods that delegate to properties() for backward compatibility

    /// Get fluid density [kg/m³] at given conditions (backward compatibility)
    fn density(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.density)
    }

    /// Get dynamic viscosity [Pa·s] at given conditions (backward compatibility)
    fn dynamic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.dynamic_viscosity)
    }

    /// Get kinematic viscosity [m²/s] at given conditions (backward compatibility)
    fn kinematic_viscosity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties(temperature, pressure)?
            .kinematic_viscosity()
    }

    /// Get specific heat capacity [J/(kg·K)] at given conditions (backward compatibility)
    fn specific_heat(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.specific_heat)
    }

    /// Get thermal conductivity [W/(m·K)] at given conditions (backward compatibility)
    fn thermal_conductivity(&self, temperature: T, pressure: T) -> Result<T, Error> {
        Ok(self.properties(temperature, pressure)?.thermal_conductivity)
    }

    /// Get Prandtl number (dimensionless) at given conditions (backward compatibility)
    fn prandtl_number(&self, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties(temperature, pressure)?.prandtl_number()
    }

    /// Calculate Reynolds number for given flow conditions (backward compatibility)
    fn reynolds_number(
        &self,
        velocity: T,
        characteristic_length: T,
        temperature: T,
        pressure: T,
    ) -> Result<T, Error> {
        self.properties(temperature, pressure)?
            .reynolds_number(velocity, characteristic_length)
    }
}

/// Type alias for simplified naming when using constant property fluids
/// This provides a shorter name for the most common fluid model in simulations
pub type Fluid<T> = ConstantPropertyFluid<T>;

/// Constant property fluid model (incompressible, Newtonian)
///
/// This model assumes fluid properties are independent of temperature and pressure.
/// Suitable for isothermal, incompressible flow simulations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ConstantPropertyFluid<T: RealField + Copy> {
    /// Descriptive name of the fluid
    pub name: String,
    /// Constant density [kg/m³]
    pub density: T,
    /// Constant dynamic viscosity [Pa·s]
    pub viscosity: T,
    /// Constant specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Constant thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
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
    /// Values from validated physics constants
    pub fn water_20c() -> Result<Self, Error>
    where
        T: FromPrimitive,
    {
        use crate::constants::physics_validated::{fluid_dynamics, thermodynamics};

        let density = T::from_f64(fluid_dynamics::WATER_DENSITY_20C).ok_or_else(|| {
            Error::InvalidInput(format!(
                "Cannot represent water density ({} kg/m³) in target type T",
                fluid_dynamics::WATER_DENSITY_20C
            ))
        })?;
        let viscosity =
            T::from_f64(fluid_dynamics::WATER_DYNAMIC_VISCOSITY_20C).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent water viscosity ({} Pa·s) in target type T",
                    fluid_dynamics::WATER_DYNAMIC_VISCOSITY_20C
                ))
            })?;
        let specific_heat =
            T::from_f64(thermodynamics::WATER_SPECIFIC_HEAT_20C).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent water specific heat ({} J/(kg·K)) in target type T",
                    thermodynamics::WATER_SPECIFIC_HEAT_20C
                ))
            })?;
        let thermal_conductivity = T::from_f64(thermodynamics::WATER_THERMAL_CONDUCTIVITY_20C)
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent water thermal conductivity ({} W/(m·K)) in target type T",
                    thermodynamics::WATER_THERMAL_CONDUCTIVITY_20C
                ))
            })?;

        Ok(Self {
            name: "Water at 20°C, 1 atm".to_string(),
            density,
            viscosity,
            specific_heat,
            thermal_conductivity,
        })
    }

    /// Create air at 20°C and 1 atm
    /// Values from validated physics constants
    pub fn air_20c() -> Result<Self, Error>
    where
        T: FromPrimitive,
    {
        use crate::constants::physics_validated::{fluid_dynamics, thermodynamics};

        let density = T::from_f64(fluid_dynamics::AIR_DENSITY_20C).ok_or_else(|| {
            Error::InvalidInput(format!(
                "Cannot represent air density ({} kg/m³) in target type T",
                fluid_dynamics::AIR_DENSITY_20C
            ))
        })?;
        let viscosity =
            T::from_f64(fluid_dynamics::AIR_DYNAMIC_VISCOSITY_20C).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent air viscosity ({} Pa·s) in target type T",
                    fluid_dynamics::AIR_DYNAMIC_VISCOSITY_20C
                ))
            })?;
        let specific_heat =
            T::from_f64(thermodynamics::AIR_SPECIFIC_HEAT_CP_20C).ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent air specific heat ({} J/(kg·K)) in target type T",
                    thermodynamics::AIR_SPECIFIC_HEAT_CP_20C
                ))
            })?;
        let thermal_conductivity = T::from_f64(thermodynamics::AIR_THERMAL_CONDUCTIVITY_20C)
            .ok_or_else(|| {
                Error::InvalidInput(format!(
                    "Cannot represent air thermal conductivity ({} W/(m·K)) in target type T",
                    thermodynamics::AIR_THERMAL_CONDUCTIVITY_20C
                ))
            })?;

        Ok(Self {
            name: "Air at 20°C, 1 atm".to_string(),
            density,
            viscosity,
            specific_heat,
            thermal_conductivity,
        })
    }

    /// Get the base constant density value
    pub fn const_density(&self) -> T {
        self.density
    }

    /// Get the base constant dynamic viscosity value
    pub fn const_dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Get the base constant kinematic viscosity
    pub fn const_kinematic_viscosity(&self) -> T {
        self.viscosity / self.density
    }

    /// Get density (alias for const_density)
    pub fn density(&self) -> T {
        self.density
    }

    /// Get viscosity (alias for const_dynamic_viscosity)
    pub fn viscosity(&self) -> T {
        self.viscosity
    }

    /// Get kinematic viscosity (alias for const_kinematic_viscosity)
    pub fn kinematic_viscosity(&self) -> T {
        self.viscosity / self.density
    }

    /// Get specific heat
    pub fn specific_heat(&self) -> T {
        self.specific_heat
    }

    /// Get thermal conductivity
    pub fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    /// Get dynamic viscosity (compatibility method)
    pub fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    /// Calculate Reynolds number for given flow conditions
    pub fn reynolds_number(&self, velocity: T, characteristic_length: T) -> T {
        self.density * velocity * characteristic_length / self.viscosity
    }
}

impl<T: RealField + Copy> crate::domains::material_properties::traits::FluidProperties<T>
    for ConstantPropertyFluid<T>
{
    fn density(&self) -> T {
        self.density
    }

    fn dynamic_viscosity(&self) -> T {
        self.viscosity
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }
}

impl<T: RealField + Copy> FluidModel<T> for ConstantPropertyFluid<T> {
    fn properties(&self, _temperature: T, _pressure: T) -> Result<FluidProperties<T>, Error> {
        Ok(FluidProperties {
            density: self.density,
            dynamic_viscosity: self.viscosity,
            specific_heat: self.specific_heat,
            thermal_conductivity: self.thermal_conductivity,
        })
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
    fn properties(&self, temperature: T, pressure: T) -> Result<FluidProperties<T>, Error> {
        if temperature <= T::zero() {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }

        // Calculate density using ideal gas law
        let r_universal = T::from_f64(Self::R_UNIVERSAL).ok_or_else(|| {
            Error::InvalidInput(
                "Cannot represent universal gas constant in target type T".to_string(),
            )
        })?;
        let r_specific = r_universal / self.molar_mass;
        let density = pressure / (r_specific * temperature);

        // Calculate viscosity using Sutherland's formula
        let t_ratio = temperature / self.temperature_ref;
        let t_ratio_3_2 = t_ratio.powf(T::from_f64(1.5).ok_or_else(|| {
            Error::InvalidInput("Cannot represent 1.5 in target type T".to_string())
        })?);
        let numerator = self.temperature_ref + self.sutherland_constant;
        let denominator = temperature + self.sutherland_constant;
        let dynamic_viscosity = self.viscosity_ref * t_ratio_3_2 * numerator / denominator;

        Ok(FluidProperties {
            density,
            dynamic_viscosity,
            specific_heat: self.specific_heat_cp,
            thermal_conductivity: self.thermal_conductivity_coeff * temperature.sqrt(),
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

// Deprecated Fluid struct removed - use ConstantPropertyFluid with FluidModel trait
