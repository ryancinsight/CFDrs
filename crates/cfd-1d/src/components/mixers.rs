//! Mixer components for microfluidic networks

use super::{Component, constants};
use cfd_core::{Error, Result, Fluid};
use nalgebra::RealField;
use num_traits::{FromPrimitive, Float};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Micromixer component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micromixer<T: RealField> {
    /// Number of inlets
    pub n_inlets: usize,
    /// Mixing efficiency (0-1)
    pub efficiency: T,
    /// Hydraulic resistance
    pub resistance: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + FromPrimitive + Float> Micromixer<T> {
    /// Create a new micromixer
    pub fn new(n_inlets: usize, resistance: T) -> Self {
        Self {
            n_inlets,
            efficiency: T::from_f64(constants::DEFAULT_MIXING_EFFICIENCY).unwrap_or_else(T::one),
            resistance,
            parameters: HashMap::new(),
        }
    }
}

impl<T: RealField + FromPrimitive + Float> Component<T> for Micromixer<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        self.resistance.clone()
    }

    fn component_type(&self) -> &str {
        "Micromixer"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "efficiency" => self.efficiency = value.max(T::zero()).min(T::one()),
            "resistance" => self.resistance = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }
}