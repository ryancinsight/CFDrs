//! Material properties domain module
//!
//! Provides fluid and material property models following Domain-Driven Design

use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Named constants for material properties
const SOLID_LIKE_VISCOSITY: f64 = 1e6;  // High viscosity for zero shear rate
const YIELD_STRESS_VISCOSITY: f64 = 1e10; // Very high viscosity below yield stress
const DEFAULT_WATER_DENSITY: f64 = 998.2; // kg/m³ at 20°C
const DEFAULT_WATER_VISCOSITY: f64 = 1.002e-3; // Pa·s at 20°C
const DEFAULT_AIR_DENSITY: f64 = 1.225; // kg/m³ at 15°C, sea level
const DEFAULT_AIR_VISCOSITY: f64 = 1.81e-5; // Pa·s at 15°C

/// Fluid properties abstraction
pub trait FluidProperties<T: RealField>: Send + Sync {
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
pub trait SolidProperties<T: RealField>: Send + Sync {
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
pub trait InterfaceProperties<T: RealField>: Send + Sync {
    /// Get surface tension
    fn surface_tension(&self) -> T;
    
    /// Get contact angle
    fn contact_angle(&self) -> T;
    
    /// Get wetting properties
    fn wetting_properties(&self) -> WettingProperties<T>;
}

/// Wetting properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WettingProperties<T: RealField> {
    /// Contact angle with solid surface
    pub contact_angle: T,
    /// Surface energy
    pub surface_energy: T,
    /// Adhesion energy
    pub adhesion_energy: T,
}

/// Newtonian fluid implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NewtonianFluid<T: RealField> {
    /// Fluid density
    pub density: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
}

impl<T: RealField> FluidProperties<T> for NewtonianFluid<T> {
    fn density(&self) -> T {
        self.density.clone()
    }
    
    fn dynamic_viscosity(&self) -> T {
        self.viscosity.clone()
    }
    
    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity.clone()
    }
    
    fn specific_heat(&self) -> T {
        self.specific_heat.clone()
    }
}

/// Non-Newtonian fluid models
pub mod non_newtonian {
    use super::*;
    
    /// Power-law fluid model
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct PowerLawFluid<T: RealField> {
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
    
    impl<T: RealField> FluidProperties<T> for PowerLawFluid<T> {
        fn density(&self) -> T {
            self.density.clone()
        }
        
        fn dynamic_viscosity(&self) -> T {
            // For power-law fluids, we return the consistency index as base viscosity
            // Actual viscosity depends on shear rate: μ = K * γ^(n-1)
            // This should be calculated with the actual shear rate when available
            self.consistency_index.clone()
        }
        
        fn thermal_conductivity(&self) -> T {
            self.thermal_conductivity.clone()
        }
        
        fn specific_heat(&self) -> T {
            self.specific_heat.clone()
        }
    }
    
    // Extension methods for non-Newtonian fluids
    impl<T: RealField> PowerLawFluid<T> {
        pub fn dynamic_viscosity_at_shear_rate(&self, shear_rate: T) -> T {
            // Power-law model: μ = K * γ^(n-1)
            // where K is consistency index, n is flow behavior index, γ is shear rate
            if shear_rate > T::zero() {
                self.consistency_index.clone() * shear_rate.powf(self.flow_behavior_index.clone() - T::one())
            } else {
                // At zero shear rate, use a large viscosity to represent solid-like behavior
                self.consistency_index.clone() * T::from_f64(SOLID_LIKE_VISCOSITY).unwrap_or(T::from_f64(SOLID_LIKE_VISCOSITY).unwrap())
            }
        }
    }
    
    /// Bingham plastic fluid model
    #[derive(Debug, Clone, Serialize, Deserialize)]
    pub struct BinghamFluid<T: RealField> {
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
    
    impl<T: RealField> FluidProperties<T> for BinghamFluid<T> {
        fn density(&self) -> T {
            self.density.clone()
        }
        
        fn dynamic_viscosity(&self) -> T {
            // For Bingham plastics, return plastic viscosity as base
            // Actual behavior depends on shear stress vs yield stress
            self.plastic_viscosity.clone()
        }
        
        fn thermal_conductivity(&self) -> T {
            self.thermal_conductivity.clone()
        }
        
        fn specific_heat(&self) -> T {
            self.specific_heat.clone()
        }
    }
    
    // Extension methods for Bingham fluids
    impl<T: RealField> BinghamFluid<T> {
        pub fn dynamic_viscosity_at_shear_stress(&self, shear_stress: T) -> T {
            // Bingham model: 
            // If τ < τ_y: material behaves as solid (infinite viscosity)
            // If τ ≥ τ_y: μ = μ_p + τ_y/γ
            let shear_stress_abs = shear_stress.abs();
            if shear_stress_abs < self.yield_stress {
                // Below yield stress - solid-like behavior
                T::from_f64(YIELD_STRESS_VISCOSITY).unwrap_or(T::from_f64(YIELD_STRESS_VISCOSITY).unwrap())
            } else {
                // Above yield stress - flows with plastic viscosity
                // Effective viscosity includes yield stress contribution
                // μ_eff = μ_p + τ_y/γ where γ = (τ - τ_y)/μ_p
                let shear_rate = (shear_stress_abs - self.yield_stress.clone()) / self.plastic_viscosity.clone();
                if shear_rate > T::zero() {
                    self.plastic_viscosity.clone() + self.yield_stress.clone() / shear_rate
                } else {
                    self.plastic_viscosity.clone()
                }
            }
        }
    }
}

/// Elastic solid implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ElasticSolid<T: RealField> {
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

impl<T: RealField> SolidProperties<T> for ElasticSolid<T> {
    fn density(&self) -> T {
        self.density.clone()
    }
    
    fn youngs_modulus(&self) -> T {
        self.youngs_modulus.clone()
    }
    
    fn poissons_ratio(&self) -> T {
        self.poissons_ratio.clone()
    }
    
    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity.clone()
    }
    
    fn specific_heat(&self) -> T {
        self.specific_heat.clone()
    }
    
    fn thermal_expansion(&self) -> T {
        self.thermal_expansion.clone()
    }
}

/// Fluid-solid interface implementation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FluidSolidInterface<T: RealField> {
    /// Surface tension
    pub surface_tension: T,
    /// Contact angle
    pub contact_angle: T,
    /// Wetting properties
    pub wetting: WettingProperties<T>,
}

impl<T: RealField> InterfaceProperties<T> for FluidSolidInterface<T> {
    fn surface_tension(&self) -> T {
        self.surface_tension.clone()
    }
    
    fn contact_angle(&self) -> T {
        self.contact_angle.clone()
    }
    
    fn wetting_properties(&self) -> WettingProperties<T> {
        self.wetting.clone()
    }
}

/// Material property database
pub struct MaterialDatabase<T: RealField> {
    /// Fluid properties database
    fluids: HashMap<String, Box<dyn FluidProperties<T>>>,
    /// Solid properties database
    solids: HashMap<String, Box<dyn SolidProperties<T>>>,
    /// Interface properties database
    interfaces: HashMap<String, Box<dyn InterfaceProperties<T>>>,
}

impl<T: RealField> MaterialDatabase<T> {
    /// Create new material database
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
    pub fn get_fluid(&self, name: &str) -> Option<&dyn FluidProperties<T>> {
        self.fluids.get(name).map(|f| f.as_ref())
    }
    
    /// Get solid properties by name
    pub fn get_solid(&self, name: &str) -> Option<&dyn SolidProperties<T>> {
        self.solids.get(name).map(|s| s.as_ref())
    }
    
    /// Get interface properties by name
    pub fn get_interface(&self, name: &str) -> Option<&dyn InterfaceProperties<T>> {
        self.interfaces.get(name).map(|i| i.as_ref())
    }
    
    /// List available fluids
    pub fn list_fluids(&self) -> Vec<&str> {
        self.fluids.keys().map(|s| s.as_str()).collect()
    }
    
    /// List available solids
    pub fn list_solids(&self) -> Vec<&str> {
        self.solids.keys().map(|s| s.as_str()).collect()
    }
    
    /// List available interfaces
    pub fn list_interfaces(&self) -> Vec<&str> {
        self.interfaces.keys().map(|s| s.as_str()).collect()
    }
}

/// Material properties service following Domain Service pattern
pub struct MaterialPropertiesService<T: RealField> {
    /// Material database
    database: MaterialDatabase<T>,
    /// Property calculators
    calculators: HashMap<String, Box<dyn PropertyCalculator<T>>>,
}

/// Property calculator abstraction
pub trait PropertyCalculator<T: RealField>: Send + Sync {
    /// Calculate derived property
    fn calculate(&self, properties: &HashMap<String, T>) -> Result<T, String>;
    
    /// Get calculator name
    fn name(&self) -> &str;
    
    /// Get required input properties
    fn required_properties(&self) -> Vec<&str>;
}

impl<T: RealField> MaterialPropertiesService<T> {
    /// Create new material properties service
    pub fn new() -> Self {
        Self {
            database: MaterialDatabase::new(),
            calculators: HashMap::new(),
        }
    }
    
    /// Get material database
    pub fn database(&self) -> &MaterialDatabase<T> {
        &self.database
    }
    
    /// Get mutable material database
    pub fn database_mut(&mut self) -> &mut MaterialDatabase<T> {
        &mut self.database
    }
    
    /// Register property calculator
    pub fn register_calculator(&mut self, name: String, calculator: Box<dyn PropertyCalculator<T>>) {
        self.calculators.insert(name, calculator);
    }
    
    /// Calculate derived property
    pub fn calculate_property(&self, calculator_name: &str, properties: &HashMap<String, T>) -> Result<T, String> {
        if let Some(calculator) = self.calculators.get(calculator_name) {
            calculator.calculate(properties)
        } else {
            Err(format!("Calculator '{}' not found", calculator_name))
        }
    }
}

impl<T: RealField> Default for MaterialDatabase<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField> Default for MaterialPropertiesService<T> {
    fn default() -> Self {
        Self::new()
    }
}
