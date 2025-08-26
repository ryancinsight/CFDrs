//! Factory for creating microfluidic components

use super::{
    constants, CircularChannel, Component, HashMap, Micropump, Microvalve, RectangularChannel,
};
use cfd_core::{Error, Result};
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
            _ => Err(Error::InvalidConfiguration(format!(
                "Unknown component type: {component_type}"
            ))),
        }
    }
}
