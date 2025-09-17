//! Material properties service for managing material data

use super::database::MaterialDatabase;
use super::property_calculators::PropertyCalculator;
use nalgebra::RealField;
use std::collections::HashMap;

/// Material properties service
pub struct MaterialPropertiesService<T: RealField + Copy> {
    /// Material database
    pub database: MaterialDatabase<T>,
    /// Property calculators
    pub calculators: HashMap<String, Box<dyn PropertyCalculator<T>>>,
}

impl<T: RealField + Copy> MaterialPropertiesService<T> {
    /// Create new material properties service
    #[must_use]
    pub fn new() -> Self {
        let mut service = Self {
            database: MaterialDatabase::new(),
            calculators: HashMap::new(),
        };

        // Register default calculators
        use super::property_calculators::{
            KinematicViscosityCalculator, PrandtlNumberCalculator, ReynoldsNumberCalculator,
        };

        service.register_calculator(
            "kinematic_viscosity".to_string(),
            Box::new(KinematicViscosityCalculator),
        );
        service.register_calculator(
            "reynolds_number".to_string(),
            Box::new(ReynoldsNumberCalculator),
        );
        service.register_calculator(
            "prandtl_number".to_string(),
            Box::new(PrandtlNumberCalculator),
        );

        service
    }

    /// Register property calculator
    pub fn register_calculator(
        &mut self,
        name: String,
        calculator: Box<dyn PropertyCalculator<T>>,
    ) {
        self.calculators.insert(name, calculator);
    }

    /// Calculate property using named calculator
    #[must_use] pub fn calculate_property(
        &self,
        calculator_name: &str,
        properties: &HashMap<String, T>,
    ) -> Option<T> {
        self.calculators
            .get(calculator_name)
            .and_then(|calc| calc.calculate(properties))
    }

    /// Get calculator by name
    #[must_use]
    pub fn get_calculator(&self, name: &str) -> Option<&dyn PropertyCalculator<T>> {
        self.calculators.get(name).map(std::convert::AsRef::as_ref)
    }

    /// List available calculators
    #[must_use]
    pub fn list_calculators(&self) -> Vec<String> {
        self.calculators.keys().cloned().collect()
    }
}

impl<T: RealField + Copy> Default for MaterialPropertiesService<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_service_initialization() {
        let service = MaterialPropertiesService::<f64>::new();

        // Check that default calculators are registered
        assert!(service.get_calculator("kinematic_viscosity").is_some());
        assert!(service.get_calculator("reynolds_number").is_some());
        assert!(service.get_calculator("prandtl_number").is_some());
    }

    #[test]
    fn test_property_calculation() {
        let service = MaterialPropertiesService::<f64>::new();

        let mut properties = HashMap::new();
        properties.insert("density".to_string(), 1000.0);
        properties.insert("dynamic_viscosity".to_string(), 0.001);

        let kinematic_viscosity = service
            .calculate_property("kinematic_viscosity", &properties)
            .unwrap();
        assert_eq!(kinematic_viscosity, 1e-6);
    }
}
