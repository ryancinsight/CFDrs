//! Material properties domain module
//!
//! Provides fluid and material property models following Domain-Driven Design

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Named constants for material properties
const SOLID_LIKE_VISCOSITY: f64 = 1e6; // High viscosity for zero shear rate
const YIELD_STRESS_VISCOSITY: f64 = 1e10; // Very high viscosity below yield stress

/// Fluid properties abstraction
pub trait FluidProperties<T: RealField + Copy>: Send + Sync {
    /// Get density
    fn density(&self) -> T;

    /// Get dynamic viscosity
    fn dynamic_viscosity(&self) -> T;

    /// Get kinematic viscosity
    fn kinematic_viscosity(&self) -> T {
        self.dynamic_viscosity() / self.density()
    }

    /// Get thermal conductivity
    fn thermal_conductivity(&self) -> T;

    /// Get specific heat capacity
    fn specific_heat(&self) -> T;

    /// Get thermal diffusivity
    fn thermal_diffusivity(&self) -> T {
        self.thermal_conductivity() / (self.density() * self.specific_heat())
    }

    /// Get Prandtl number
    fn prandtl_number(&self) -> T {
        self.kinematic_viscosity() / self.thermal_diffusivity()
    }
}

/// Solid properties abstraction
pub trait SolidProperties<T: RealField + Copy>: Send + Sync {
    /// Get density
    fn density(&self) -> T;

    /// Get Young's modulus
    fn youngs_modulus(&self) -> T;

    /// Get Poisson's ratio
    fn poissons_ratio(&self) -> T;

    /// Get thermal conductivity
    fn thermal_conductivity(&self) -> T;

    /// Get specific heat capacity
    fn specific_heat(&self) -> T;

    /// Get thermal expansion coefficient
    fn thermal_expansion(&self) -> T;
}

/// Interface properties abstraction
pub trait InterfaceProperties<T: RealField + Copy>: Send + Sync {
    /// Get surface tension
    fn surface_tension(&self) -> T;

    /// Get contact angle
    fn contact_angle(&self) -> T;

    /// Get wetting properties
    fn wetting_properties(&self) -> WettingProperties<T>;
}

/// Wetting properties
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct WettingProperties<T: RealField + Copy> {
    /// Contact angle with solid surface
    pub contact_angle: T,
    /// Surface energy
    pub surface_energy: T,
    /// Adhesion energy
    pub adhesion_energy: T,
}

/// Currenttonian fluid implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CurrenttonianFluid<T: RealField + Copy> {
    /// Fluid density
    pub density: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
}

impl<T: RealField + Copy> FluidProperties<T> for CurrenttonianFluid<T> {
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

/// Non-Currenttonian fluid models
pub mod non_newtonian {
    use super::{
        Deserialize, FluidProperties, RealField, Serialize, SOLID_LIKE_VISCOSITY,
        YIELD_STRESS_VISCOSITY,
    };

    /// Power-law fluid model
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct PowerLawFluid<T: RealField + Copy> {
        /// Fluid density
        pub density: T,
        /// Consistency index
        pub consistency_index: T,
        /// Flow behavior index
        pub flow_behavior_index: T,
        /// Thermal conductivity
        pub thermal_conductivity: T,
        /// Specific heat capacity
        pub specific_heat: T,
    }

    impl<T: RealField + Copy> FluidProperties<T> for PowerLawFluid<T> {
        fn density(&self) -> T {
            self.density
        }

        fn dynamic_viscosity(&self) -> T {
            // For power-law fluids, we return the consistency index as base viscosity
            // Actual viscosity depends on shear rate: μ = K * γ^(n-1)
            // This should be calculated with the actual shear rate when available
            self.consistency_index
        }

        fn thermal_conductivity(&self) -> T {
            self.thermal_conductivity
        }

        fn specific_heat(&self) -> T {
            self.specific_heat
        }
    }

    // Extension methods for non-Currenttonian fluids
    impl<T: RealField + Copy> PowerLawFluid<T> {
        pub fn dynamic_viscosity_at_shear_rate(&self, shear_rate: T) -> T {
            // Power-law model: μ = K * γ^(n-1)
            // where K is consistency index, n is flow behavior index, γ is shear rate
            if shear_rate > T::zero() {
                self.consistency_index * shear_rate.powf(self.flow_behavior_index - T::one())
            } else {
                // At zero shear rate, use a large viscosity to represent solid-like behavior
                self.consistency_index
                    * T::from_f64(SOLID_LIKE_VISCOSITY)
                        .unwrap_or(T::from_f64(SOLID_LIKE_VISCOSITY).unwrap_or_else(|| T::one()))
            }
        }
    }

    /// Bingham plastic fluid model
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct BinghamFluid<T: RealField + Copy> {
        /// Fluid density
        pub density: T,
        /// Plastic viscosity
        pub plastic_viscosity: T,
        /// Yield stress
        pub yield_stress: T,
        /// Thermal conductivity
        pub thermal_conductivity: T,
        /// Specific heat capacity
        pub specific_heat: T,
    }

    impl<T: RealField + Copy> FluidProperties<T> for BinghamFluid<T> {
        fn density(&self) -> T {
            self.density
        }

        fn dynamic_viscosity(&self) -> T {
            // For Bingham plastics, return plastic viscosity as base
            // Actual behavior depends on shear stress vs yield stress
            self.plastic_viscosity
        }

        fn thermal_conductivity(&self) -> T {
            self.thermal_conductivity
        }

        fn specific_heat(&self) -> T {
            self.specific_heat
        }
    }

    // Extension methods for Bingham fluids
    impl<T: RealField + Copy> BinghamFluid<T> {
        pub fn dynamic_viscosity_at_shear_stress(&self, shear_stress: T) -> T {
            // Bingham model:
            // If τ < τ_y: material behaves as solid (infinite viscosity)
            // If τ ≥ τ_y: μ = μ_p + τ_y/γ
            let shear_stress_abs = shear_stress.abs();
            if shear_stress_abs < self.yield_stress {
                // Below yield stress - solid-like behavior
                T::from_f64(YIELD_STRESS_VISCOSITY)
                    .unwrap_or(T::from_f64(YIELD_STRESS_VISCOSITY).unwrap_or_else(|| T::one()))
            } else {
                // Above yield stress - flows with plastic viscosity
                // Effective viscosity includes yield stress contribution
                // μ_eff = μ_p + τ_y/γ where γ = (τ - τ_y)/μ_p
                let shear_rate = (shear_stress_abs - self.yield_stress) / self.plastic_viscosity;
                if shear_rate > T::zero() {
                    self.plastic_viscosity + self.yield_stress / shear_rate
                } else {
                    self.plastic_viscosity
                }
            }
        }
    }
}

/// Elastic solid implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElasticSolid<T: RealField + Copy> {
    /// Solid density
    pub density: T,
    /// Young's modulus
    pub youngs_modulus: T,
    /// Poisson's ratio
    pub poissons_ratio: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
    /// Thermal expansion coefficient
    pub thermal_expansion: T,
}

impl<T: RealField + Copy> SolidProperties<T> for ElasticSolid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn youngs_modulus(&self) -> T {
        self.youngs_modulus
    }

    fn poissons_ratio(&self) -> T {
        self.poissons_ratio
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }

    fn thermal_expansion(&self) -> T {
        self.thermal_expansion
    }
}

/// Fluid-solid interface implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidSolidInterface<T: RealField + Copy> {
    /// Surface tension
    pub surface_tension: T,
    /// Contact angle
    pub contact_angle: T,
    /// Wetting properties
    pub wetting: WettingProperties<T>,
}

impl<T: RealField + Copy> InterfaceProperties<T> for FluidSolidInterface<T> {
    fn surface_tension(&self) -> T {
        self.surface_tension
    }

    fn contact_angle(&self) -> T {
        self.contact_angle
    }

    fn wetting_properties(&self) -> WettingProperties<T> {
        self.wetting
    }
}

/// Material property database
pub struct MaterialDatabase<T: RealField + Copy> {
    /// Fluid properties database
    fluids: HashMap<String, Box<dyn FluidProperties<T>>>,
    /// Solid properties database
    solids: HashMap<String, Box<dyn SolidProperties<T>>>,
    /// Interface properties database
    interfaces: HashMap<String, Box<dyn InterfaceProperties<T>>>,
}

impl<T: RealField + Copy> MaterialDatabase<T> {
    /// Create new material database
    #[must_use]
    pub fn new() -> Self {
        Self {
            fluids: HashMap::new(),
            solids: HashMap::new(),
            interfaces: HashMap::new(),
        }
    }

    /// Add fluid to database
    pub fn add_fluid(&mut self, name: String, fluid: Box<dyn FluidProperties<T>>) {
        self.fluids.insert(name, fluid);
    }

    /// Add solid to database
    pub fn add_solid(&mut self, name: String, solid: Box<dyn SolidProperties<T>>) {
        self.solids.insert(name, solid);
    }

    /// Add interface to database
    pub fn add_interface(&mut self, name: String, interface: Box<dyn InterfaceProperties<T>>) {
        self.interfaces.insert(name, interface);
    }

    /// Get fluid properties by name
    #[must_use]
    pub fn get_fluid(&self, name: &str) -> Option<&dyn FluidProperties<T>> {
        self.fluids.get(name).map(std::convert::AsRef::as_ref)
    }

    /// Get solid properties by name
    #[must_use]
    pub fn get_solid(&self, name: &str) -> Option<&dyn SolidProperties<T>> {
        self.solids.get(name).map(std::convert::AsRef::as_ref)
    }

    /// Get interface properties by name
    #[must_use]
    pub fn get_interface(&self, name: &str) -> Option<&dyn InterfaceProperties<T>> {
        self.interfaces.get(name).map(std::convert::AsRef::as_ref)
    }

    /// List available fluids
    #[must_use]
    pub fn list_fluids(&self) -> Vec<&str> {
        self.fluids
            .keys()
            .map(std::string::String::as_str)
            .collect()
    }

    /// List available solids
    #[must_use]
    pub fn list_solids(&self) -> Vec<&str> {
        self.solids
            .keys()
            .map(std::string::String::as_str)
            .collect()
    }

    /// List available interfaces
    #[must_use]
    pub fn list_interfaces(&self) -> Vec<&str> {
        self.interfaces
            .keys()
            .map(std::string::String::as_str)
            .collect()
    }
}

/// Material properties service following Domain Service pattern
pub struct MaterialPropertiesService<T: RealField + Copy> {
    /// Material database
    database: MaterialDatabase<T>,
    /// Property calculators
    calculators: HashMap<String, Box<dyn PropertyCalculator<T>>>,
}

/// Property calculator abstraction
pub trait PropertyCalculator<T: RealField + Copy>: Send + Sync {
    /// Calculate derived property
    fn calculate(&self, properties: &HashMap<String, T>) -> Result<T, String>;

    /// Get calculator name
    fn name(&self) -> &str;

    /// Get required input properties
    fn required_properties(&self) -> Vec<&str>;
}

/// Kinematic viscosity calculator: ν = μ/ρ
pub struct KinematicViscosityCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for KinematicViscosityCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Result<T, String> {
        let viscosity = properties
            .get("dynamic_viscosity")
            .ok_or_else(|| "Missing dynamic_viscosity".to_string())?;
        let density = properties
            .get("density")
            .ok_or_else(|| "Missing density".to_string())?;

        if *density <= T::zero() {
            return Err("Density must be positive".to_string());
        }

        Ok(*viscosity / *density)
    }

    fn name(&self) -> &str {
        "kinematic_viscosity"
    }

    fn required_properties(&self) -> Vec<&str> {
        vec!["dynamic_viscosity", "density"]
    }
}

/// Reynolds number calculator: Re = ρVL/μ
pub struct ReynoldsNumberCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for ReynoldsNumberCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Result<T, String> {
        let density = properties
            .get("density")
            .ok_or_else(|| "Missing density".to_string())?;
        let velocity = properties
            .get("velocity")
            .ok_or_else(|| "Missing velocity".to_string())?;
        let length = properties
            .get("characteristic_length")
            .ok_or_else(|| "Missing characteristic_length".to_string())?;
        let viscosity = properties
            .get("dynamic_viscosity")
            .ok_or_else(|| "Missing dynamic_viscosity".to_string())?;

        if *viscosity <= T::zero() {
            return Err("Viscosity must be positive".to_string());
        }

        Ok(*density * *velocity * *length / *viscosity)
    }

    fn name(&self) -> &str {
        "reynolds_number"
    }

    fn required_properties(&self) -> Vec<&str> {
        vec![
            "density",
            "velocity",
            "characteristic_length",
            "dynamic_viscosity",
        ]
    }
}

/// Prandtl number calculator: Pr = μCp/k
pub struct PrandtlNumberCalculator;

impl<T: RealField + Copy> PropertyCalculator<T> for PrandtlNumberCalculator {
    fn calculate(&self, properties: &HashMap<String, T>) -> Result<T, String> {
        let viscosity = properties
            .get("dynamic_viscosity")
            .ok_or_else(|| "Missing dynamic_viscosity".to_string())?;
        let cp = properties
            .get("specific_heat_cp")
            .ok_or_else(|| "Missing specific_heat_cp".to_string())?;
        let conductivity = properties
            .get("thermal_conductivity")
            .ok_or_else(|| "Missing thermal_conductivity".to_string())?;

        if *conductivity <= T::zero() {
            return Err("Thermal conductivity must be positive".to_string());
        }

        Ok(*viscosity * *cp / *conductivity)
    }

    fn name(&self) -> &str {
        "prandtl_number"
    }

    fn required_properties(&self) -> Vec<&str> {
        vec![
            "dynamic_viscosity",
            "specific_heat_cp",
            "thermal_conductivity",
        ]
    }
}

impl<T: RealField + Copy> MaterialPropertiesService<T> {
    /// Create new material properties service
    #[must_use]
    pub fn new() -> Self {
        let mut service = Self {
            database: MaterialDatabase::new(),
            calculators: HashMap::new(),
        };

        // Register standard calculators
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

    /// Get material database
    #[must_use]
    pub fn database(&self) -> &MaterialDatabase<T> {
        &self.database
    }

    /// Get mutable material database
    pub fn database_mut(&mut self) -> &mut MaterialDatabase<T> {
        &mut self.database
    }

    /// Register property calculator
    pub fn register_calculator(
        &mut self,
        name: String,
        calculator: Box<dyn PropertyCalculator<T>>,
    ) {
        self.calculators.insert(name, calculator);
    }

    /// Calculate derived property
    pub fn calculate_property(
        &self,
        calculator_name: &str,
        properties: &HashMap<String, T>,
    ) -> Result<T, String> {
        if let Some(calculator) = self.calculators.get(calculator_name) {
            calculator.calculate(properties)
        } else {
            Err(format!("Calculator '{calculator_name}' not found"))
        }
    }
}

impl<T: RealField + Copy> Default for MaterialDatabase<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> Default for MaterialPropertiesService<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {}
