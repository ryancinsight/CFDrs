//! Basic fluid properties and calculations

use crate::error::Error;
use crate::physics::fluid::thermophysical;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, KinematicViscosity, Length, MassDensity, SpecificHeatCapacity,
    ThermalConductivity, ThermalDiffusivity, ThermodynamicTemperature, Velocity,
};
use eunomia::NumericElement;
use eunomia::RealField;
use serde::{Deserialize, Serialize};

/// A block of computed fluid properties at a single state point
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct FluidProperties<T> {
    /// Density [kg/m³]
    pub density: MassDensity<T>,
    /// Dynamic viscosity [Pa·s]
    pub dynamic_viscosity: DynamicViscosity<T>,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: SpecificHeatCapacity<T>,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: ThermalConductivity<T>,
}

impl<T: RealField + Copy> FluidProperties<T> {
    /// Create new fluid properties
    pub fn new(
        density: MassDensity<T>,
        dynamic_viscosity: DynamicViscosity<T>,
        specific_heat: SpecificHeatCapacity<T>,
        thermal_conductivity: ThermalConductivity<T>,
    ) -> Self {
        Self {
            density,
            dynamic_viscosity,
            specific_heat,
            thermal_conductivity,
        }
    }
}

impl<T: RealField + NumericElement + Copy> FluidProperties<T> {
    /// Calculate kinematic viscosity [m²/s] from base properties
    ///
    /// # Errors
    /// Returns an error if density is non-positive
    pub fn kinematic_viscosity(&self) -> Result<KinematicViscosity<T>, Error> {
        if self.density.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput("Density must be positive".to_string()));
        }
        Ok(self.dynamic_viscosity / self.density)
    }

    /// Calculate Prandtl number from base properties
    ///
    /// # Errors
    /// Returns an error if thermal conductivity is non-positive
    pub fn prandtl_number(&self) -> Result<Dimensionless<T>, Error> {
        if self.thermal_conductivity.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Thermal conductivity must be positive".to_string(),
            ));
        }
        Ok(self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity)
    }

    /// Calculate Reynolds number for given flow conditions
    ///
    /// # Errors
    /// Returns an error if the dynamic viscosity is zero or negative
    pub fn reynolds_number(
        &self,
        velocity: Velocity<T>,
        characteristic_length: Length<T>,
    ) -> Result<Dimensionless<T>, Error> {
        if self.dynamic_viscosity.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Viscosity must be positive".to_string(),
            ));
        }
        Ok(self.density * velocity * characteristic_length / self.dynamic_viscosity)
    }

    /// Calculate Peclet number for given flow conditions
    ///
    /// # Errors
    /// Returns an error if viscosity, thermal conductivity, density, or specific heat are invalid
    pub fn peclet_number(
        &self,
        velocity: Velocity<T>,
        characteristic_length: Length<T>,
    ) -> Result<Dimensionless<T>, Error> {
        let re = self.reynolds_number(velocity, characteristic_length)?;
        let pr = self.prandtl_number()?;
        Ok(re * pr)
    }

    /// Calculate thermal diffusivity [m²/s]
    ///
    /// # Errors
    /// Returns an error if density, specific heat, or thermal conductivity
    /// violates the Proteus thermophysical-property contract.
    pub fn thermal_diffusivity(&self) -> Result<ThermalDiffusivity<T>, Error> {
        thermophysical::thermal_diffusivity(
            self.density,
            self.specific_heat,
            self.thermal_conductivity,
        )
    }

    /// Calculate speed of sound for ideal gas \[m/s]
    /// Requires ratio of specific heats (gamma) as parameter
    ///
    /// # Errors
    /// Returns an error if any input parameter is non-positive
    pub fn speed_of_sound(
        &self,
        gamma: Dimensionless<T>,
        temperature: ThermodynamicTemperature<T>,
        gas_constant: SpecificHeatCapacity<T>,
    ) -> Result<Velocity<T>, Error> {
        if temperature.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Temperature must be positive".to_string(),
            ));
        }
        Ok(Velocity::from_base(<T as NumericElement>::sqrt(
            gamma.into_base() * gas_constant.into_base() * temperature.into_base(),
        )))
    }

    /// Calculate Mach number for given flow velocity
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidInput` if sound speed is non-positive.
    pub fn mach_number(
        &self,
        velocity: Velocity<T>,
        sound_speed: Velocity<T>,
    ) -> Result<Dimensionless<T>, Error> {
        if sound_speed.into_base() <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidInput(
                "Sound speed must be positive".to_string(),
            ));
        }
        Ok(Dimensionless::from_base(
            velocity.into_base() / sound_speed.into_base(),
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn water() -> FluidProperties<f64> {
        FluidProperties::new(
            MassDensity::from_base(1_000.0),
            DynamicViscosity::from_base(0.001),
            SpecificHeatCapacity::from_base(4_186.0),
            ThermalConductivity::from_base(0.6),
        )
    }

    #[test]
    fn typed_metrics_match_independent_formula_values() {
        let properties = water();
        let velocity = Velocity::from_base(2.0);
        let length = Length::from_base(0.01);

        let kinematic = properties
            .kinematic_viscosity()
            .expect("positive density produces kinematic viscosity")
            .into_base();
        let prandtl = properties
            .prandtl_number()
            .expect("positive conductivity produces Prandtl number")
            .into_base();
        let reynolds = properties
            .reynolds_number(velocity, length)
            .expect("positive viscosity produces Reynolds number")
            .into_base();
        let peclet = properties
            .peclet_number(velocity, length)
            .expect("positive thermophysical properties produce Peclet number")
            .into_base();
        let diffusivity = properties
            .thermal_diffusivity()
            .expect("positive thermophysical properties produce diffusivity")
            .into_base();
        let speed_of_sound = properties
            .speed_of_sound(
                Dimensionless::from_base(1.4),
                ThermodynamicTemperature::from_base(300.0),
                SpecificHeatCapacity::from_base(287.05),
            )
            .expect("positive ideal-gas inputs produce sound speed")
            .into_base();
        let mach = properties
            .mach_number(Velocity::from_base(2.0), Velocity::from_base(343.2))
            .expect("positive sound speed produces Mach number")
            .into_base();

        let expected_kinematic = 0.001 / 1_000.0;
        let expected_prandtl = 0.001 * 4_186.0 / 0.6;
        let expected_reynolds = 1_000.0 * 2.0 * 0.01 / 0.001;
        let expected_peclet = expected_reynolds * expected_prandtl;
        let expected_diffusivity = 0.6 / (1_000.0 * 4_186.0);
        let expected_speed_of_sound = (1.4_f64 * 287.05 * 300.0).sqrt();
        let expected_mach = 2.0 / 343.2;

        // Each path has at most four native f64 roundings; the bound is the
        // independently evaluated result scaled by that error budget.
        for (actual, expected) in [
            (kinematic, expected_kinematic),
            (prandtl, expected_prandtl),
            (reynolds, expected_reynolds),
            (peclet, expected_peclet),
            (diffusivity, expected_diffusivity),
            (speed_of_sound, expected_speed_of_sound),
            (mach, expected_mach),
        ] {
            let tolerance = 4.0 * f64::EPSILON * expected.abs().max(1.0);
            assert!((actual - expected).abs() <= tolerance);
        }
    }

    #[test]
    fn typed_metric_validation_reports_the_invalid_quantity() {
        let mut properties = water();
        properties.density = MassDensity::from_base(-1.0);
        let error = properties
            .kinematic_viscosity()
            .expect_err("negative density is outside the metric domain");
        assert_eq!(error.to_string(), "Invalid input: Density must be positive");

        let mut properties = water();
        properties.thermal_conductivity = ThermalConductivity::from_base(0.0);
        let error = properties
            .prandtl_number()
            .expect_err("zero conductivity is outside the metric domain");
        assert_eq!(
            error.to_string(),
            "Invalid input: Thermal conductivity must be positive"
        );
    }
}
