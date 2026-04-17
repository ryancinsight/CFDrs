//! Factory for creating microfluidic components

use super::{
    constants, try_real_from_f64, CircularChannel, Component, HashMap, Micropump, Microvalve,
    OrganCompartment, PorousMembrane, RectangularChannel,
};
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Factory for creating microfluidic components
pub struct ComponentFactory;

impl ComponentFactory {
    fn required_param<T: RealField + Copy + FromPrimitive + Float>(
        params: &HashMap<String, T>,
        name: &str,
    ) -> Result<T> {
        params
            .get(name)
            .copied()
            .ok_or_else(|| Error::InvalidConfiguration(format!("Missing {name} parameter")))
    }

    fn optional_param<T: RealField + Copy + FromPrimitive + Float>(
        params: &HashMap<String, T>,
        name: &str,
        default: f64,
        label: &str,
    ) -> Result<T> {
        if let Some(&v) = params.get(name) {
            Ok(v)
        } else {
            try_real_from_f64(default, label)
        }
    }

    /// Create a component from type string and parameters
    pub fn create<T: RealField + Copy + FromPrimitive + Float>(
        component_type: &str,
        params: &HashMap<String, T>,
    ) -> Result<Box<dyn Component<T>>> {
        match component_type {
            "RectangularChannel" => {
                let length = Self::required_param(params, "length")?;
                let width = Self::required_param(params, "width")?;
                let height = Self::required_param(params, "height")?;
                let roughness = Self::optional_param(
                    params,
                    "roughness",
                    constants::DEFAULT_ROUGHNESS,
                    "default roughness",
                )?;
                Ok(Box::new(RectangularChannel::new(
                    length, width, height, roughness,
                )))
            }
            "CircularChannel" => {
                let length = Self::required_param(params, "length")?;
                let diameter = Self::required_param(params, "diameter")?;
                let roughness = Self::optional_param(
                    params,
                    "roughness",
                    constants::DEFAULT_ROUGHNESS,
                    "default roughness",
                )?;
                Ok(Box::new(CircularChannel::new(length, diameter, roughness)))
            }
            "Micropump" => {
                let max_flow_rate = Self::required_param(params, "max_flow_rate")?;
                let max_pressure = Self::required_param(params, "max_pressure")?;
                Ok(Box::new(Micropump::new(max_flow_rate, max_pressure)))
            }
            "Microvalve" => {
                let cv = Self::optional_param(
                    params,
                    "cv",
                    constants::DEFAULT_VALVE_CV,
                    "default valve cv",
                )?;
                Ok(Box::new(Microvalve::new(cv)))
            }
            "PorousMembrane" => {
                let thickness = Self::required_param(params, "thickness")?;
                let width = Self::required_param(params, "width")?;
                let height = Self::required_param(params, "height")?;
                let pore_radius = Self::required_param(params, "pore_radius")?;
                let porosity = Self::required_param(params, "porosity")?;
                Ok(Box::new(PorousMembrane::new(
                    thickness,
                    width,
                    height,
                    pore_radius,
                    porosity,
                )))
            }
            "OrganCompartment" => {
                let length = Self::required_param(params, "length")?;
                let width = Self::required_param(params, "width")?;
                let height = Self::required_param(params, "height")?;
                let hydraulic_resistance = Self::required_param(params, "hydraulic_resistance")?;
                Ok(Box::new(OrganCompartment::new(
                    length,
                    width,
                    height,
                    hydraulic_resistance,
                )))
            }
            _ => Err(Error::InvalidConfiguration(format!(
                "Unknown component type: {component_type}"
            ))),
        }
    }
}
