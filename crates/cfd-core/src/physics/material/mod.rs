//! Material management system for multi-physics simulations
//!
//! This module provides a unified system for managing fluid, solid, and interface
//! properties, following Domain-Driven Design principles with deep vertical integration.

pub mod interface;
pub mod solid;
pub mod traits;

pub use interface::{FluidSolidInterface, WettingProperties};
pub use solid::ElasticSolid;
pub use traits::{InterfaceProperties, SolidProperties};

use crate::physics::fluid::traits::Fluid;
use nalgebra::RealField;
use std::collections::HashMap;

/// Unified material database for multi-physics simulations
pub struct MaterialDatabase<T: RealField + Copy> {
    /// Fluid materials (using the unified Fluid trait)
    fluids: HashMap<String, Box<dyn Fluid<T>>>,
    /// Solid materials
    solids: HashMap<String, Box<dyn SolidProperties<T>>>,
    /// Interface properties
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

    /// Add fluid material
    pub fn add_fluid(&mut self, name: String, fluid: Box<dyn Fluid<T>>) {
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
    pub fn get_fluid(&self, name: &str) -> Option<&dyn Fluid<T>> {
        self.fluids.get(name).map(|f| f.as_ref())
    }

    /// Get solid material by name
    #[must_use]
    pub fn get_solid(&self, name: &str) -> Option<&dyn SolidProperties<T>> {
        self.solids.get(name).map(|s| s.as_ref())
    }

    /// Get interface properties by name
    #[must_use]
    pub fn get_interface(&self, name: &str) -> Option<&dyn InterfaceProperties<T>> {
        self.interfaces.get(name).map(|i| i.as_ref())
    }

    /// Initialize with common materials
    #[must_use]
    pub fn with_common_materials() -> Self
    where
        T: nalgebra::RealField + Copy + num_traits::FromPrimitive,
    {
        let mut db = Self::new();

        // Add common fluids from the physics::fluid database
        if let Ok(water) = crate::physics::fluid::database::water_20c::<T>() {
            db.add_fluid("water".to_string(), Box::new(water));
        }
        if let Ok(air) = crate::physics::fluid::database::air_20c::<T>() {
            db.add_fluid("air".to_string(), Box::new(air));
        }

        // Add common solids
        db.add_solid("steel".to_string(), Box::new(ElasticSolid::steel()));
        db.add_solid("aluminum".to_string(), Box::new(ElasticSolid::aluminum()));

        // Add common interfaces
        db.add_interface(
            "water_air".to_string(),
            Box::new(FluidSolidInterface::water_air()),
        );

        db
    }
}

impl<T: RealField + Copy> Default for MaterialDatabase<T> {
    fn default() -> Self {
        Self::new()
    }
}
