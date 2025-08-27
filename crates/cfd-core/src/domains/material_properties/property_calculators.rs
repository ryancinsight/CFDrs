//! Property calculators for derived material properties

use nalgebra::RealField;
use std::collections::HashMap;

/// Property calculator trait
pub trait PropertyCalculator<T: RealField + Copy>: Send + Sync {
    /// Calculate property value
    fn calculate(&self, properties: &HashMap<String, T>) -> Option<T>;

    /// Get required property names
    fn required_properties(&self) -> Vec<String>;

    /// Get calculator name
    fn name(&self) -> &str;
}

/// Kinematic viscosity calculator
pub struct KinematicViscosityCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for KinematicViscosityCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Option<T> {
        let density = properties.get("density")?;
        let dynamic_viscosity = properties.get("dynamic_viscosity")?;

        if *density > T::zero() {
            Some(*dynamic_viscosity / *density)
        } else {
            None
        }
    }

    fn required_properties(&self) -> Vec<String> {
        vec!["density".to_string(), "dynamic_viscosity".to_string()]
    }

    fn name(&self) -> &str {
        "Kinematic Viscosity"
    }
}

/// Reynolds number calculator
pub struct ReynoldsNumberCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for ReynoldsNumberCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Option<T> {
        let density = properties.get("density")?;
        let velocity = properties.get("velocity")?;
        let length_scale = properties.get("length_scale")?;
        let dynamic_viscosity = properties.get("dynamic_viscosity")?;

        if *dynamic_viscosity > T::zero() {
            Some(*density * *velocity * *length_scale / *dynamic_viscosity)
        } else {
            None
        }
    }

    fn required_properties(&self) -> Vec<String> {
        vec![
            "density".to_string(),
            "velocity".to_string(),
            "length_scale".to_string(),
            "dynamic_viscosity".to_string(),
        ]
    }

    fn name(&self) -> &str {
        "Reynolds Number"
    }
}

/// Prandtl number calculator
pub struct PrandtlNumberCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for PrandtlNumberCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Option<T> {
        let specific_heat = properties.get("specific_heat")?;
        let dynamic_viscosity = properties.get("dynamic_viscosity")?;
        let thermal_conductivity = properties.get("thermal_conductivity")?;

        if *thermal_conductivity > T::zero() {
            Some(*specific_heat * *dynamic_viscosity / *thermal_conductivity)
        } else {
            None
        }
    }

    fn required_properties(&self) -> Vec<String> {
        vec![
            "specific_heat".to_string(),
            "dynamic_viscosity".to_string(),
            "thermal_conductivity".to_string(),
        ]
    }

    fn name(&self) -> &str {
        "Prandtl Number"
    }
}
