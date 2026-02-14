//! Membrane and organ components for organ-on-chip style 1D simulations.

use super::Component;
use crate::resistance::models::{FlowConditions, MembranePoreModel, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::physics::fluid::Fluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Porous membrane represented as many cylindrical pores in parallel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PorousMembrane<T: RealField + Copy> {
    /// Membrane thickness [m]
    pub thickness: T,
    /// Membrane width [m]
    pub width: T,
    /// Membrane height [m]
    pub height: T,
    /// Pore radius [m]
    pub pore_radius: T,
    /// Fraction of open pore area (0..1)
    pub porosity: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> PorousMembrane<T> {
    /// Create a new porous membrane component.
    pub fn new(thickness: T, width: T, height: T, pore_radius: T, porosity: T) -> Self {
        Self {
            thickness,
            width,
            height,
            pore_radius,
            porosity,
            parameters: HashMap::new(),
        }
    }

    /// Membrane frontal area [m²]
    pub fn area(&self) -> T {
        self.width * self.height
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for PorousMembrane<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        let model = MembranePoreModel::new(
            self.thickness,
            self.width,
            self.height,
            self.pore_radius,
            self.porosity,
        );
        let conditions = FlowConditions::new(T::zero());
        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| T::from_f64(1e12).unwrap_or_else(T::zero))
    }

    fn component_type(&self) -> &'static str {
        "PorousMembrane"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "thickness" => self.thickness = value,
            "width" => self.width = value,
            "height" => self.height = value,
            "pore_radius" => self.pore_radius = value,
            "porosity" => self.porosity = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.thickness * self.area())
    }
}

/// Organ chamber represented as a compartment with user-defined hydraulic resistance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganCompartment<T: RealField + Copy> {
    /// Chamber length [m]
    pub length: T,
    /// Chamber width [m]
    pub width: T,
    /// Chamber height [m]
    pub height: T,
    /// Lumped hydraulic resistance [Pa·s/m³]
    pub hydraulic_resistance: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> OrganCompartment<T> {
    /// Create a new organ compartment component.
    pub fn new(length: T, width: T, height: T, hydraulic_resistance: T) -> Self {
        Self {
            length,
            width,
            height,
            hydraulic_resistance,
            parameters: HashMap::new(),
        }
    }

    /// Chamber area [m²]
    pub fn area(&self) -> T {
        self.width * self.height
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for OrganCompartment<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        self.hydraulic_resistance
    }

    fn component_type(&self) -> &'static str {
        "OrganCompartment"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "width" => self.width = value,
            "height" => self.height = value,
            "hydraulic_resistance" => self.hydraulic_resistance = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.length * self.area())
    }
}
