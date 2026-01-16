//! Validation utilities for fluid properties
//!
//! Ensures physical consistency and reasonable bounds for all fluid properties

use super::FluidProperties;
use crate::error::Error;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Physical bounds for fluid properties
pub struct PropertyBounds<T: RealField + Copy> {
    /// Minimum density [kg/m³]
    pub density_min: T,
    /// Maximum density [kg/m³]
    pub density_max: T,
    /// Minimum viscosity [Pa·s]
    pub viscosity_min: T,
    /// Maximum viscosity [Pa·s]
    pub viscosity_max: T,
    /// Minimum specific heat [J/(kg·K)]
    pub specific_heat_min: T,
    /// Maximum specific heat [J/(kg·K)]
    pub specific_heat_max: T,
    /// Minimum thermal conductivity [W/(m·K)]
    pub thermal_conductivity_min: T,
    /// Maximum thermal conductivity [W/(m·K)]
    pub thermal_conductivity_max: T,
}

impl<T: RealField + FromPrimitive + Copy> Default for PropertyBounds<T> {
    fn default() -> Self {
        Self {
            // Covers gases to heavy liquids
            density_min: T::from_f64(0.01).unwrap_or_else(T::zero),
            density_max: T::from_f64(20000.0).unwrap_or_else(T::one),
            // From superfluid helium to highly viscous materials
            viscosity_min: T::from_f64(1e-7).unwrap_or_else(T::zero),
            viscosity_max: T::from_f64(1e6).unwrap_or_else(T::one),
            // Reasonable range for most fluids
            specific_heat_min: T::from_f64(100.0).unwrap_or_else(T::zero),
            specific_heat_max: T::from_f64(10000.0).unwrap_or_else(T::one),
            // From insulators to liquid metals
            thermal_conductivity_min: T::from_f64(0.001).unwrap_or_else(T::zero),
            thermal_conductivity_max: T::from_f64(1000.0).unwrap_or_else(T::one),
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
            "Density {} is outside valid range [{}, {}]",
            properties.density, bounds.density_min, bounds.density_max
        )));
    }

    // Check viscosity
    if properties.dynamic_viscosity < bounds.viscosity_min
        || properties.dynamic_viscosity > bounds.viscosity_max
    {
        return Err(Error::InvalidInput(format!(
            "Viscosity {} is outside valid range [{}, {}]",
            properties.dynamic_viscosity, bounds.viscosity_min, bounds.viscosity_max
        )));
    }

    // Check specific heat
    if properties.specific_heat < bounds.specific_heat_min
        || properties.specific_heat > bounds.specific_heat_max
    {
        return Err(Error::InvalidInput(format!(
            "Specific heat {} is outside valid range [{}, {}]",
            properties.specific_heat, bounds.specific_heat_min, bounds.specific_heat_max
        )));
    }

    // Check thermal conductivity
    if properties.thermal_conductivity < bounds.thermal_conductivity_min
        || properties.thermal_conductivity > bounds.thermal_conductivity_max
    {
        return Err(Error::InvalidInput(format!(
            "Thermal conductivity {} is outside valid range [{}, {}]",
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
pub fn validate_reynolds<T: RealField + Copy>(reynolds: T) -> Result<(), Error> {
    if reynolds < T::zero() {
        return Err(Error::InvalidInput(
            "Reynolds number cannot be negative".to_string(),
        ));
    }

    // Warn for extreme Reynolds numbers
    let re_max = T::from_f64(1e8).unwrap_or_else(T::one);
    if reynolds > re_max {
        return Err(Error::InvalidInput(format!(
            "Reynolds number {reynolds} exceeds typical range (may indicate error)"
        )));
    }

    Ok(())
}

/// Check Prandtl number validity
///
/// # Errors
///
/// Returns `Error::InvalidInput` if Prandtl number is non-positive or outside typical range (0.001-100000).
pub fn validate_prandtl<T: RealField + Copy>(prandtl: T) -> Result<(), Error> {
    if prandtl <= T::zero() {
        return Err(Error::InvalidInput(
            "Prandtl number must be positive".to_string(),
        ));
    }

    // Typical range: 0.001 (liquid metals) to 100000 (heavy oils)
    let pr_min = T::from_f64(0.001).unwrap_or_else(T::zero);
    let pr_max = T::from_f64(100_000.0).unwrap_or_else(T::one);

    if prandtl < pr_min || prandtl > pr_max {
        return Err(Error::InvalidInput(format!(
            "Prandtl number {prandtl} is outside typical range [{pr_min}, {pr_max}]"
        )));
    }

    Ok(())
}

/// Validate temperature for physical reasonableness
///
/// # Errors
///
/// TODO: Returns `Error::InvalidInput` if temperature is non-positive (at/below absolute zero) or exceeds 10000 K.
/// DEPENDENCIES: Implement context-aware temperature validation based on application
/// BLOCKED BY: Fixed validation limits don't account for different physics regimes
/// PRIORITY: Medium - Important for specialized applications (plasma, cryogenics)
/// Returns `Error::InvalidInput` if temperature is non-positive (at/below absolute zero) or exceeds 10000 K.
pub fn validate_temperature<T: RealField + Copy>(temperature: T) -> Result<(), Error> {
    let t_min = T::zero(); // Absolute zero
    let t_max = T::from_f64(10000.0).unwrap_or_else(T::one); // Reasonable upper limit

    if temperature <= t_min {
        return Err(Error::InvalidInput(
            "Temperature must be positive (above absolute zero)".to_string(),
        ));
    }

    if temperature > t_max {
        return Err(Error::InvalidInput(format!(
            "Temperature {temperature} K exceeds reasonable range"
        )));
    }

    Ok(())
}

/// Validate pressure for physical reasonableness
///
/// # Errors
///
/// Returns `Error::InvalidInput` if pressure is non-positive or exceeds 1 GPa.
pub fn validate_pressure<T: RealField + Copy>(pressure: T) -> Result<(), Error> {
    if pressure <= T::zero() {
        return Err(Error::InvalidInput("Pressure must be positive".to_string()));
    }

    // Maximum pressure: ~1 GPa (deep ocean/industrial)
    let p_max = T::from_f64(1e9).unwrap_or_else(T::one);
    if pressure > p_max {
        return Err(Error::InvalidInput(format!(
            "Pressure {pressure} Pa exceeds typical range"
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
            1000.0, // density
            0.001,  // viscosity
            4186.0, // specific heat
            0.6,    // thermal conductivity
        );
        assert!(validate_properties(&valid, &bounds).is_ok());

        // Invalid density
        let invalid_density = FluidProperties::new(
            -100.0, // negative density
            0.001, 4186.0, 0.6,
        );
        assert!(validate_properties(&invalid_density, &bounds).is_err());
    }

    #[test]
    fn test_reynolds_validation() {
        assert!(validate_reynolds(1000.0).is_ok());
        assert!(validate_reynolds(0.0).is_ok());
        assert!(validate_reynolds(-100.0).is_err());
        assert!(validate_reynolds(1e9).is_err()); // Too large
    }
}
