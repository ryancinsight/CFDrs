//! Validation utilities for fluid properties
//!
//! Ensures physical consistency and reasonable bounds for all fluid properties

use super::FluidProperties;
use crate::error::Error;
use aequitas::systems::si::quantities::{
    DynamicViscosity, MassDensity, SpecificHeatCapacity, ThermalConductivity,
};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};

/// Physical bounds for fluid properties
pub struct PropertyBounds<T> {
    /// Minimum density [kg/m³]
    pub density_min: MassDensity<T>,
    /// Maximum density [kg/m³]
    pub density_max: MassDensity<T>,
    /// Minimum viscosity [Pa·s]
    pub viscosity_min: DynamicViscosity<T>,
    /// Maximum viscosity [Pa·s]
    pub viscosity_max: DynamicViscosity<T>,
    /// Minimum specific heat [J/(kg·K)]
    pub specific_heat_min: SpecificHeatCapacity<T>,
    /// Maximum specific heat [J/(kg·K)]
    pub specific_heat_max: SpecificHeatCapacity<T>,
    /// Minimum thermal conductivity [W/(m·K)]
    pub thermal_conductivity_min: ThermalConductivity<T>,
    /// Maximum thermal conductivity [W/(m·K)]
    pub thermal_conductivity_max: ThermalConductivity<T>,
}

impl<T: RealField + FloatElement + Copy> Default for PropertyBounds<T> {
    fn default() -> Self {
        Self {
            // Covers gases to heavy liquids
            density_min: MassDensity::from_base(<T as FloatElement>::from_f64(0.01)),
            density_max: MassDensity::from_base(<T as FloatElement>::from_f64(20000.0)),
            // From superfluid helium to highly viscous materials
            viscosity_min: DynamicViscosity::from_base(<T as FloatElement>::from_f64(1e-7)),
            viscosity_max: DynamicViscosity::from_base(<T as FloatElement>::from_f64(1e6)),
            // Reasonable range for most fluids
            specific_heat_min: SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(
                100.0,
            )),
            specific_heat_max: SpecificHeatCapacity::from_base(<T as FloatElement>::from_f64(
                10000.0,
            )),
            // From insulators to liquid metals
            thermal_conductivity_min: ThermalConductivity::from_base(
                <T as FloatElement>::from_f64(0.001),
            ),
            thermal_conductivity_max: ThermalConductivity::from_base(
                <T as FloatElement>::from_f64(1000.0),
            ),
        }
    }
}

/// Validate fluid properties against physical bounds
///
/// # Errors
///
/// Returns `Error::InvalidInput` if any property is outside the specified bounds.
pub fn validate_properties<T: RealField + Copy>(
    properties: &FluidProperties<T>,
    bounds: &PropertyBounds<T>,
) -> Result<(), Error> {
    // Check density
    if properties.density < bounds.density_min || properties.density > bounds.density_max {
        return Err(Error::InvalidInput(format!(
            "Density {:?} is outside valid range [{:?}, {:?}]",
            properties.density, bounds.density_min, bounds.density_max
        )));
    }

    // Check viscosity
    if properties.dynamic_viscosity < bounds.viscosity_min
        || properties.dynamic_viscosity > bounds.viscosity_max
    {
        return Err(Error::InvalidInput(format!(
            "Viscosity {:?} is outside valid range [{:?}, {:?}]",
            properties.dynamic_viscosity, bounds.viscosity_min, bounds.viscosity_max
        )));
    }

    // Check specific heat
    if properties.specific_heat < bounds.specific_heat_min
        || properties.specific_heat > bounds.specific_heat_max
    {
        return Err(Error::InvalidInput(format!(
            "Specific heat {:?} is outside valid range [{:?}, {:?}]",
            properties.specific_heat, bounds.specific_heat_min, bounds.specific_heat_max
        )));
    }

    // Check thermal conductivity
    if properties.thermal_conductivity < bounds.thermal_conductivity_min
        || properties.thermal_conductivity > bounds.thermal_conductivity_max
    {
        return Err(Error::InvalidInput(format!(
            "Thermal conductivity {:?} is outside valid range [{:?}, {:?}]",
            properties.thermal_conductivity,
            bounds.thermal_conductivity_min,
            bounds.thermal_conductivity_max
        )));
    }

    Ok(())
}

/// Check dimensionless number validity
///
/// # Errors
///
/// Returns `Error::InvalidInput` if Reynolds number is negative or exceeds typical range (>1e8).
pub fn validate_reynolds<T: RealField + FloatElement + Copy>(reynolds: T) -> Result<(), Error> {
    if reynolds < <T as NumericElement>::ZERO {
        return Err(Error::InvalidInput(
            "Reynolds number cannot be negative".to_string(),
        ));
    }

    // Warn for extreme Reynolds numbers
    let re_max = <T as FloatElement>::from_f64(1e8);
    if reynolds > re_max {
        return Err(Error::InvalidInput(format!(
            "Reynolds number {reynolds:?} exceeds typical range (may indicate error)"
        )));
    }

    Ok(())
}

/// Check Prandtl number validity
///
/// # Errors
///
/// Returns `Error::InvalidInput` if Prandtl number is non-positive or outside typical range (0.001-100000).
pub fn validate_prandtl<T: RealField + FloatElement + Copy>(prandtl: T) -> Result<(), Error> {
    if prandtl <= <T as NumericElement>::ZERO {
        return Err(Error::InvalidInput(
            "Prandtl number must be positive".to_string(),
        ));
    }

    // Typical range: 0.001 (liquid metals) to 100000 (heavy oils)
    let pr_min = <T as FloatElement>::from_f64(0.001);
    let pr_max = <T as FloatElement>::from_f64(100_000.0);

    if prandtl < pr_min || prandtl > pr_max {
        return Err(Error::InvalidInput(format!(
            "Prandtl number {prandtl:?} is outside typical range [{pr_min:?}, {pr_max:?}]"
        )));
    }

    Ok(())
}

/// Validate temperature for physical reasonableness
///
/// # Errors
///
/// Returns `Error::InvalidInput` if temperature is non-positive (at/below absolute zero) or exceeds 10000 K.
pub fn validate_temperature<T: RealField + FloatElement + Copy>(
    temperature: T,
) -> Result<(), Error> {
    let t_min = <T as NumericElement>::ZERO; // Absolute zero
    let t_max = <T as FloatElement>::from_f64(10000.0); // Reasonable upper limit

    if temperature <= t_min {
        return Err(Error::InvalidInput(
            "Temperature must be positive (above absolute zero)".to_string(),
        ));
    }

    if temperature > t_max {
        return Err(Error::InvalidInput(format!(
            "Temperature {temperature:?} K exceeds reasonable range"
        )));
    }

    Ok(())
}

/// Validate pressure for physical reasonableness
///
/// # Errors
///
/// Returns `Error::InvalidInput` if pressure is non-positive or exceeds 1 GPa.
pub fn validate_pressure<T: RealField + FloatElement + Copy>(pressure: T) -> Result<(), Error> {
    if pressure <= <T as NumericElement>::ZERO {
        return Err(Error::InvalidInput("Pressure must be positive".to_string()));
    }

    // Maximum pressure: ~1 GPa (deep ocean/industrial)
    let p_max = <T as FloatElement>::from_f64(1e9);
    if pressure > p_max {
        return Err(Error::InvalidInput(format!(
            "Pressure {pressure:?} Pa exceeds typical range"
        )));
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_property_validation() {
        let bounds = PropertyBounds::<f64>::default();

        // Valid properties (water-like)
        let valid = FluidProperties::new(
            MassDensity::from_base(1000.0),
            DynamicViscosity::from_base(0.001),
            SpecificHeatCapacity::from_base(4186.0),
            ThermalConductivity::from_base(0.6),
        );
        validate_properties(&valid, &bounds)
            .expect("invariant: water-like fluid properties are valid");

        // Invalid density
        let invalid_density = FluidProperties::new(
            MassDensity::from_base(-100.0),
            DynamicViscosity::from_base(0.001),
            SpecificHeatCapacity::from_base(4186.0),
            ThermalConductivity::from_base(0.6),
        );
        validate_properties(&invalid_density, &bounds)
            .expect_err("invariant: negative density is rejected");
    }

    #[test]
    fn test_reynolds_validation() {
        validate_reynolds(1000.0).expect("invariant: nominal Reynolds number is valid");
        validate_reynolds(0.0).expect("invariant: zero Reynolds number is valid");
        validate_reynolds(-100.0).expect_err("invariant: negative Reynolds number is rejected");
        validate_reynolds(1e9).expect_err("invariant: excessive Reynolds number is rejected");
    }
}
