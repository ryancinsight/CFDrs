//! Sensor components for microfluidic networks

use super::{Component, constants};
use cfd_core::{Error, Result, Fluid};
use nalgebra::RealField;
use num_traits::{FromPrimitive, Float};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Flow sensor component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowSensor<T: RealField> {
    /// Sensor resistance [Pa·s/m³]
    pub resistance: T,
    /// Measurement range [m³/s]
    pub range: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + FromPrimitive + Float> FlowSensor<T> {
    /// Create a new flow sensor
    pub fn new(resistance: T, range: T) -> Self {
        Self {
            resistance,
            range,
            parameters: HashMap::new(),
        }
    }
}

impl<T: RealField + FromPrimitive + Float> Component<T> for FlowSensor<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        self.resistance.clone()
    }

    fn component_type(&self) -> &str {
        "FlowSensor"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "resistance" => self.resistance = value,
            "range" => self.range = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }
}