//! Membrane and organ components for organ-on-chip style 1D simulations.

use super::{real_from_f64, Component};
use crate::physics::resistance::models::{FlowConditions, MembranePoreModel, ResistanceModel};
use aequitas::systems::si::quantities::{Area, Length, Volume};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_core::CfdScalar;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Porous membrane represented as many cylindrical pores in parallel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PorousMembrane<T: CfdScalar + Copy> {
    /// Membrane thickness \[m]
    pub thickness: Length<T>,
    /// Membrane width \[m]
    pub width: Length<T>,
    /// Membrane height \[m]
    pub height: Length<T>,
    /// Pore radius \[m]
    pub pore_radius: Length<T>,
    /// Fraction of open pore area (0..1)
    pub porosity: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: CfdScalar + Copy + SafeFromF64> PorousMembrane<T> {
    /// Create a new porous membrane component.
    pub fn new(
        thickness: Length<T>,
        width: Length<T>,
        height: Length<T>,
        pore_radius: Length<T>,
        porosity: T,
    ) -> Self {
        Self {
            thickness,
            width,
            height,
            pore_radius,
            porosity,
            parameters: HashMap::new(),
        }
    }

    /// Membrane frontal area \[m²]
    pub fn area(&self) -> Area<T> {
        Area::from_base(self.width.into_base() * self.height.into_base())
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> Component<T> for PorousMembrane<T> {
    fn resistance(&self, fluid: &ConstantPropertyFluid<T>) -> T {
        let model = MembranePoreModel::new(
            self.thickness.into_base(),
            self.width.into_base(),
            self.height.into_base(),
            self.pore_radius.into_base(),
            self.porosity,
        );
        let conditions = FlowConditions::new(T::ZERO);
        model
            .calculate_resistance(fluid, &conditions)
            .unwrap_or_else(|_| real_from_f64(1e12))
    }

    fn component_type(&self) -> &'static str {
        "PorousMembrane"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "thickness" => self.thickness = Length::from_base(value),
            "width" => self.width = Length::from_base(value),
            "height" => self.height = Length::from_base(value),
            "pore_radius" => self.pore_radius = Length::from_base(value),
            "porosity" => self.porosity = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<Volume<T>> {
        Some(Volume::from_base(
            self.thickness.into_base() * self.area().into_base(),
        ))
    }
}

/// Organ chamber represented as a compartment with user-defined hydraulic resistance.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganCompartment<T: CfdScalar + Copy> {
    /// Chamber length \[m]
    pub length: Length<T>,
    /// Chamber width \[m]
    pub width: Length<T>,
    /// Chamber height \[m]
    pub height: Length<T>,
    /// Lumped hydraulic resistance [Pa·s/m³]
    pub hydraulic_resistance: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: CfdScalar + Copy + SafeFromF64> OrganCompartment<T> {
    /// Create a new organ compartment component.
    pub fn new(
        length: Length<T>,
        width: Length<T>,
        height: Length<T>,
        hydraulic_resistance: T,
    ) -> Self {
        Self {
            length,
            width,
            height,
            hydraulic_resistance,
            parameters: HashMap::new(),
        }
    }

    /// Chamber area \[m²]
    pub fn area(&self) -> Area<T> {
        Area::from_base(self.width.into_base() * self.height.into_base())
    }
}

impl<T: CfdScalar + Copy + SafeFromF64> Component<T> for OrganCompartment<T> {
    fn resistance(&self, _fluid: &ConstantPropertyFluid<T>) -> T {
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
            "length" => self.length = Length::from_base(value),
            "width" => self.width = Length::from_base(value),
            "height" => self.height = Length::from_base(value),
            "hydraulic_resistance" => self.hydraulic_resistance = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<Volume<T>> {
        Some(Volume::from_base(
            self.length.into_base() * self.area().into_base(),
        ))
    }
}
