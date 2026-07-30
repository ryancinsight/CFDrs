//! Newtonian fluid models with constant and variable properties

use super::thermophysical;
use super::traits::{ConstantFluid, Fluid as FluidTrait, FluidState};
use crate::error::Error;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, Pressure, SpecificHeatCapacity, TemperatureDifference,
    ThermalConductivity, ThermodynamicTemperature, Velocity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

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
    /// Constant speed of sound \[m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + Copy> ConstantPropertyFluid<T> {
    /// Create a new constant property fluid
    pub fn new(
        name: String,
        density: T,
        viscosity: T,
        specific_heat: T,
        thermal_conductivity: T,
        speed_of_sound: T,
    ) -> Self {
        Self {
            name,
            density,
            viscosity,
            specific_heat,
            thermal_conductivity,
            speed_of_sound,
        }
    }

    /// Validate that all properties are physically reasonable.
    ///
    /// Density, specific heat, and thermal conductivity are validated through the
    /// Proteus `ThermophysicalProperties` contract (finite, non-negative). Viscosity
    /// and speed of sound are checked separately since Proteus does not yet model those
    /// property kinds.
    ///
    /// # Errors
    /// Returns a descriptive error naming the first violated physical constraint.
    pub fn validate(&self) -> Result<(), Error>
    where
        T: NumericElement,
    {
        // Route the thermophysical subset through Proteus — this checks finiteness
        // and non-negativeness with descriptive error messages.
        thermophysical::validate_thermophysical_subset(
            self.density,
            self.specific_heat,
            self.thermal_conductivity,
        )?;

        if self.viscosity <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        if self.speed_of_sound <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Speed of sound must be positive".to_string(),
            ));
        }
        Ok(())
    }

    /// Water at 20°C and standard pressure
    /// Properties from NIST webbook
    ///
    /// # Errors
    /// Returns an error if validation rejects the constructed physical constants.
    pub fn water_20c() -> Result<Self, Error>
    where
        T: FloatElement,
    {
        let fluid = Self::new(
            "Water (20°C)".to_string(),
            <T as FloatElement>::from_f64(998.2), // density [kg/m³]
            <T as FloatElement>::from_f64(0.001_002), // viscosity [Pa·s]
            <T as FloatElement>::from_f64(4186.0), // specific heat [J/(kg·K)]
            <T as FloatElement>::from_f64(0.599), // thermal conductivity [W/(m·K)]
            <T as FloatElement>::from_f64(1482.0), // speed of sound [m/s]
        );
        fluid.validate()?;
        Ok(fluid)
    }

    /// Air at 20°C and standard pressure
    /// Properties from NIST webbook
    ///
    /// # Errors
    /// Returns an error if validation rejects the constructed physical constants.
    pub fn air_20c() -> Result<Self, Error>
    where
        T: FloatElement,
    {
        let fluid = Self::new(
            "Air (20°C)".to_string(),
            <T as FloatElement>::from_f64(1.204), // density [kg/m³]
            <T as FloatElement>::from_f64(1.825e-5), // viscosity [Pa·s]
            <T as FloatElement>::from_f64(1005.0), // specific heat [J/(kg·K)]
            <T as FloatElement>::from_f64(0.02538), // thermal conductivity [W/(m·K)]
            <T as FloatElement>::from_f64(343.2), // speed of sound [m/s]
        );
        fluid.validate()?;
        Ok(fluid)
    }
}

impl<T: RealField + Copy> FluidTrait<T> for ConstantPropertyFluid<T> {
    fn properties_at(&self, _temperature: T, _pressure: T) -> Result<FluidState<T>, Error> {
        Ok(FluidState {
            density: MassDensity::from_base(self.density),
            dynamic_viscosity: DynamicViscosity::from_base(self.viscosity),
            specific_heat: SpecificHeatCapacity::from_base(self.specific_heat),
            thermal_conductivity: ThermalConductivity::from_base(self.thermal_conductivity),
            speed_of_sound: Velocity::from_base(self.speed_of_sound),
        })
    }

    fn name(&self) -> &str {
        &self.name
    }
}

impl<T: RealField + Copy> ConstantFluid<T> for ConstantPropertyFluid<T> {
    fn density(&self) -> MassDensity<T> {
        MassDensity::from_base(self.density)
    }

    fn dynamic_viscosity(&self) -> DynamicViscosity<T> {
        DynamicViscosity::from_base(self.viscosity)
    }

    fn specific_heat(&self) -> SpecificHeatCapacity<T> {
        SpecificHeatCapacity::from_base(self.specific_heat)
    }

    fn thermal_conductivity(&self) -> ThermalConductivity<T> {
        ThermalConductivity::from_base(self.thermal_conductivity)
    }

    fn speed_of_sound(&self) -> Velocity<T> {
        Velocity::from_base(self.speed_of_sound)
    }
}

/// Ideal gas model with constant specific heats
///
/// Properties vary with temperature and pressure according to ideal gas law
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IdealGas<T> {
    /// Gas name
    pub name: String,
    /// Specific gas constant [J/(kg·K)]
    pub gas_constant: SpecificHeatCapacity<T>,
    /// Specific heat at constant pressure [J/(kg·K)]
    pub cp: SpecificHeatCapacity<T>,
    /// Reference viscosity [Pa·s] at reference temperature
    pub mu_ref: DynamicViscosity<T>,
    /// Reference temperature \[K] for viscosity
    pub t_ref: ThermodynamicTemperature<T>,
    /// Sutherland constant \[K] for viscosity model
    pub sutherland_constant: TemperatureDifference<T>,
    /// Thermal conductivity coefficient [J/(kg·K)]
    pub k_coeff: SpecificHeatCapacity<T>,
}

impl<T: RealField + FloatElement + Copy> IdealGas<T> {
    /// Create a new ideal gas
    pub fn new(
        name: String,
        gas_constant: SpecificHeatCapacity<T>,
        cp: SpecificHeatCapacity<T>,
        mu_ref: DynamicViscosity<T>,
        t_ref: ThermodynamicTemperature<T>,
        sutherland_constant: TemperatureDifference<T>,
    ) -> Self {
        // Thermal conductivity from Prandtl number assumption (Pr ≈ 0.7 for air)
        let pr = <T as FloatElement>::from_f64(0.7);
        let k_coeff = SpecificHeatCapacity::from_base(cp.into_base() / pr);

        Self {
            name,
            gas_constant,
            cp,
            mu_ref,
            t_ref,
            sutherland_constant,
            k_coeff,
        }
    }

    /// Calculate density from ideal gas law
    fn calculate_density(
        &self,
        temperature: ThermodynamicTemperature<T>,
        pressure: Pressure<T>,
    ) -> MassDensity<T> {
        MassDensity::from_base(
            pressure.into_base() / (self.gas_constant.into_base() * temperature.into_base()),
        )
    }

    /// Calculate viscosity using Sutherland's law
    fn calculate_viscosity(&self, temperature: ThermodynamicTemperature<T>) -> DynamicViscosity<T> {
        let temperature_base = temperature.into_base();
        let t_ref = self.t_ref.into_base();
        let sutherland_constant = self.sutherland_constant.into_base();
        let t_ratio = temperature_base / t_ref;
        let numerator = t_ref + sutherland_constant;
        let denominator = temperature_base + sutherland_constant;

        DynamicViscosity::from_base(
            self.mu_ref.into_base()
                * <T as FloatElement>::powf(t_ratio, <T as FloatElement>::from_f64(1.5))
                * numerator
                / denominator,
        )
    }

    /// Calculate thermal conductivity using kinetic theory
    fn calculate_thermal_conductivity(
        &self,
        viscosity: DynamicViscosity<T>,
    ) -> ThermalConductivity<T> {
        ThermalConductivity::from_base(viscosity.into_base() * self.k_coeff.into_base())
    }

    /// Calculate ideal-gas speed of sound from the constant-heat-capacity model.
    fn calculate_speed_of_sound(&self, temperature: ThermodynamicTemperature<T>) -> Velocity<T> {
        // c = sqrt(gamma * R * T), gamma = cp / (cp - R).
        let cp = self.cp.into_base();
        let gas_constant = self.gas_constant.into_base();
        let gamma = cp / (cp - gas_constant);
        Velocity::from_base(<T as NumericElement>::sqrt(
            gamma * gas_constant * temperature.into_base(),
        ))
    }
}

impl<T: RealField + FloatElement + Copy> FluidTrait<T> for IdealGas<T> {
    fn properties_at(&self, temperature: T, pressure: T) -> Result<FluidState<T>, Error> {
        let temperature = ThermodynamicTemperature::from_base(temperature);
        let pressure = Pressure::from_base(pressure);

        if temperature.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        if pressure.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput("Pressure must be positive".to_string()));
        }

        let density = self.calculate_density(temperature, pressure);
        let viscosity = self.calculate_viscosity(temperature);
        let thermal_conductivity = self.calculate_thermal_conductivity(viscosity);
        let speed_of_sound = self.calculate_speed_of_sound(temperature);

        Ok(FluidState {
            density,
            dynamic_viscosity: viscosity,
            specific_heat: self.cp,
            thermal_conductivity,
            speed_of_sound,
        })
    }

    fn name(&self) -> &str {
        &self.name
    }

    fn is_temperature_dependent(&self) -> bool {
        true
    }

    fn is_pressure_dependent(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base_fluid() -> ConstantPropertyFluid<f64> {
        ConstantPropertyFluid::new(
            "test".to_string(),
            998.2,     // density [kg/m³]
            0.001_002, // viscosity [Pa·s]
            4186.0,    // specific heat [J/(kg·K)]
            0.599,     // thermal conductivity [W/(m·K)]
            1482.0,    // speed of sound [m/s]
        )
    }

    #[test]
    fn water_20c_is_physically_valid() {
        ConstantPropertyFluid::<f64>::water_20c().expect("NIST water properties pass validation");
    }

    #[test]
    fn air_20c_is_physically_valid() {
        ConstantPropertyFluid::<f64>::air_20c().expect("NIST air properties pass validation");
    }

    #[test]
    fn validate_rejects_negative_density_with_proteus_message() {
        let mut fluid = base_fluid();
        fluid.density = -1.0;
        let err = fluid.validate().expect_err("negative density is invalid");
        match err {
            Error::InvalidInput(msg) => assert!(
                msg.contains("MassDensity"),
                "Proteus error must name the violated property; got: {msg}"
            ),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_nan_density_with_proteus_message() {
        let mut fluid = base_fluid();
        fluid.density = f64::NAN;
        let err = fluid.validate().expect_err("NaN density is invalid");
        match err {
            Error::InvalidInput(msg) => assert!(msg.contains("MassDensity")),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_negative_specific_heat_with_proteus_message() {
        let mut fluid = base_fluid();
        fluid.specific_heat = -1.0;
        let err = fluid
            .validate()
            .expect_err("negative specific heat is invalid");
        match err {
            Error::InvalidInput(msg) => assert!(
                msg.contains("SpecificHeatCapacity"),
                "Proteus error must name the violated property; got: {msg}"
            ),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_negative_thermal_conductivity_with_proteus_message() {
        let mut fluid = base_fluid();
        fluid.thermal_conductivity = -0.1;
        let err = fluid
            .validate()
            .expect_err("negative thermal conductivity is invalid");
        match err {
            Error::InvalidInput(msg) => assert!(
                msg.contains("ThermalConductivity"),
                "Proteus error must name the violated property; got: {msg}"
            ),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_zero_viscosity_with_descriptive_message() {
        let mut fluid = base_fluid();
        fluid.viscosity = 0.0;
        let err = fluid.validate().expect_err("zero viscosity is invalid");
        match err {
            Error::InvalidInput(msg) => {
                assert!(msg.to_lowercase().contains("viscosity"), "got: {msg}");
            }
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn validate_rejects_zero_speed_of_sound_with_descriptive_message() {
        let mut fluid = base_fluid();
        fluid.speed_of_sound = 0.0;
        let err = fluid
            .validate()
            .expect_err("zero speed of sound is invalid");
        match err {
            Error::InvalidInput(msg) => assert!(msg.to_lowercase().contains("speed"), "got: {msg}"),
            other => panic!("unexpected error variant: {other:?}"),
        }
    }

    #[test]
    fn ideal_gas_preserves_typed_state_metrics() {
        let gas = IdealGas::new(
            "Air".to_string(),
            SpecificHeatCapacity::from_base(287.0),
            SpecificHeatCapacity::from_base(1005.0),
            DynamicViscosity::from_base(1.716e-5),
            ThermodynamicTemperature::from_base(273.15),
            TemperatureDifference::from_base(110.4),
        );

        let state = <IdealGas<f64> as FluidTrait<f64>>::properties_at(&gas, 288.15, 101_325.0)
            .expect("positive ideal-gas state is valid");

        assert!((state.density.into_base() - 1.225).abs() < 0.01);
        assert!(state.dynamic_viscosity.into_base() > 0.0);
        assert!(state.thermal_conductivity.into_base() > 0.0);
        assert!((state.speed_of_sound.into_base() - 340.3).abs() < 1.0);
        assert_eq!(state.specific_heat, gas.cp);
    }

    #[test]
    fn ideal_gas_rejects_nonpositive_state_inputs() {
        let gas = IdealGas::new(
            "Air".to_string(),
            SpecificHeatCapacity::from_base(287.0),
            SpecificHeatCapacity::from_base(1005.0),
            DynamicViscosity::from_base(1.716e-5),
            ThermodynamicTemperature::from_base(273.15),
            TemperatureDifference::from_base(110.4),
        );

        assert!(<IdealGas<f64> as FluidTrait<f64>>::properties_at(&gas, 0.0, 101_325.0).is_err());
        assert!(<IdealGas<f64> as FluidTrait<f64>>::properties_at(&gas, 288.15, 0.0).is_err());
    }
}
