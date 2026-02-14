//! Factory for creating microfluidic components

use super::{
    constants, CircularChannel, Component, HashMap, Micropump, Microvalve, OrganCompartment,
    PorousMembrane, RectangularChannel,
};
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Factory for creating microfluidic components
pub struct ComponentFactory;

impl ComponentFactory {
    /// Create a component from type string and parameters
    pub fn create<T: RealField + Copy + FromPrimitive + Float>(
        component_type: &str,
        params: &HashMap<String, T>,
    ) -> Result<Box<dyn Component<T>>> {
        match component_type {
            "RectangularChannel" => {
                let length = params.get("length").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing length parameter".into())
                })?;
                let width = params
                    .get("width")
                    .ok_or_else(|| Error::InvalidConfiguration("Missing width parameter".into()))?;
                let height = params.get("height").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing height parameter".into())
                })?;
                let roughness = params.get("roughness").copied().unwrap_or_else(|| {
                    T::from_f64(constants::DEFAULT_ROUGHNESS).unwrap_or_else(|| T::zero())
                });

                Ok(Box::new(RectangularChannel::new(
                    *length, *width, *height, roughness,
                )))
            }
            "CircularChannel" => {
                let length = params.get("length").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing length parameter".into())
                })?;
                let diameter = params.get("diameter").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing diameter parameter".into())
                })?;
                let roughness = params.get("roughness").copied().unwrap_or_else(|| {
                    T::from_f64(constants::DEFAULT_ROUGHNESS).unwrap_or_else(|| T::zero())
                });

                Ok(Box::new(CircularChannel::new(
                    *length, *diameter, roughness,
                )))
            }
            "Micropump" => {
                let max_flow_rate = params.get("max_flow_rate").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing max_flow_rate parameter".into())
                })?;
                let max_pressure = params.get("max_pressure").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing max_pressure parameter".into())
                })?;

                Ok(Box::new(Micropump::new(*max_flow_rate, *max_pressure)))
            }
            "Microvalve" => {
                let cv = params.get("cv").copied().unwrap_or_else(|| {
                    T::from_f64(constants::DEFAULT_VALVE_CV).unwrap_or_else(|| T::zero())
                });

                Ok(Box::new(Microvalve::new(cv)))
            }
            "PorousMembrane" => {
                let thickness = params.get("thickness").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing thickness parameter".into())
                })?;
                let width = params
                    .get("width")
                    .ok_or_else(|| Error::InvalidConfiguration("Missing width parameter".into()))?;
                let height = params.get("height").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing height parameter".into())
                })?;
                let pore_radius = params.get("pore_radius").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing pore_radius parameter".into())
                })?;
                let porosity = params.get("porosity").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing porosity parameter".into())
                })?;

                Ok(Box::new(PorousMembrane::new(
                    *thickness,
                    *width,
                    *height,
                    *pore_radius,
                    *porosity,
                )))
            }
            "OrganCompartment" => {
                let length = params.get("length").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing length parameter".into())
                })?;
                let width = params
                    .get("width")
                    .ok_or_else(|| Error::InvalidConfiguration("Missing width parameter".into()))?;
                let height = params.get("height").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing height parameter".into())
                })?;
                let hydraulic_resistance = params.get("hydraulic_resistance").ok_or_else(|| {
                    Error::InvalidConfiguration("Missing hydraulic_resistance parameter".into())
                })?;

                Ok(Box::new(OrganCompartment::new(
                    *length,
                    *width,
                    *height,
                    *hydraulic_resistance,
                )))
            }
            _ => Err(Error::InvalidConfiguration(format!(
                "Unknown component type: {component_type}"
            ))),
        }
    }
}
