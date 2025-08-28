//! Material database for storing and retrieving material properties

use super::traits::{FluidProperties, InterfaceProperties, SolidProperties};
use nalgebra::RealField;
use std::collections::HashMap;

/// Material database for storing material properties
pub struct MaterialDatabase<T: RealField + Copy> {
    /// Fluid materials
    pub fluids: HashMap<String, Box<dyn FluidProperties<T>>>,
    /// Solid materials
    pub solids: HashMap<String, Box<dyn SolidProperties<T>>>,
    /// Interface properties
    pub interfaces: HashMap<String, Box<dyn InterfaceProperties<T>>>,
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

    /// Add fluid material
    pub fn add_fluid(&mut self, name: String, fluid: Box<dyn FluidProperties<T>>) {
        self.fluids.insert(name, fluid);
    }

    /// Add solid material
    pub fn add_solid(&mut self, name: String, solid: Box<dyn SolidProperties<T>>) {
        self.solids.insert(name, solid);
    }

    /// Add interface properties
    pub fn add_interface(&mut self, name: String, interface: Box<dyn InterfaceProperties<T>>) {
        self.interfaces.insert(name, interface);
    }

    /// Get fluid material by name
    #[must_use]
    pub fn get_fluid(&self, name: &str) -> Option<&dyn FluidProperties<T>> {
        self.fluids.get(name).map(std::convert::AsRef::as_ref)
    }

    /// Get solid material by name
    #[must_use]
    pub fn get_solid(&self, name: &str) -> Option<&dyn SolidProperties<T>> {
        self.solids.get(name).map(std::convert::AsRef::as_ref)
    }

    /// Get interface properties by name
    #[must_use]
    pub fn get_interface(&self, name: &str) -> Option<&dyn InterfaceProperties<T>> {
        self.interfaces.get(name).map(std::convert::AsRef::as_ref)
    }

    /// List available fluid materials
    #[must_use]
    pub fn list_fluids(&self) -> Vec<String> {
        self.fluids.keys().cloned().collect()
    }

    /// List available solid materials
    #[must_use]
    pub fn list_solids(&self) -> Vec<String> {
        self.solids.keys().cloned().collect()
    }

    /// List available interfaces
    #[must_use]
    pub fn list_interfaces(&self) -> Vec<String> {
        self.interfaces.keys().cloned().collect()
    }

    /// Initialize with common materials
    pub fn with_common_materials() -> Self
    where
        T: From<f64>,
    {
        let mut db = Self::new();

        // Add common fluids
        use super::fluids::NewtonianFluid;
        db.add_fluid(
            "water".to_string(),
            Box::new(NewtonianFluid::<T>::water()),
        );
        db.add_fluid("air".to_string(), Box::new(NewtonianFluid::<T>::air()));
        db.add_fluid("oil".to_string(), Box::new(NewtonianFluid::<T>::oil()));

        // Add common solids
        use super::solids::ElasticSolid;
        db.add_solid("steel".to_string(), Box::new(ElasticSolid::<T>::steel()));
        db.add_solid(
            "aluminum".to_string(),
            Box::new(ElasticSolid::<T>::aluminum()),
        );
        db.add_solid(
            "concrete".to_string(),
            Box::new(ElasticSolid::<T>::concrete()),
        );

        // Add common interfaces
        use super::interfaces::FluidSolidInterface;
        db.add_interface(
            "water_air".to_string(),
            Box::new(FluidSolidInterface::<T>::water_air()),
        );
        db.add_interface(
            "oil_water".to_string(),
            Box::new(FluidSolidInterface::<T>::oil_water()),
        );

        db
    }
}

impl<T: RealField + Copy> Default for MaterialDatabase<T> {
    fn default() -> Self {
        Self::new()
    }
}
