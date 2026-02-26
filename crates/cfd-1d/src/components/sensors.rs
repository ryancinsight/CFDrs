//! Sensor components for microfluidic networks
//!
//! # FlowSensor Resistance Theorem
//!
//! A flow sensor is modelled as an **ideal measurement device**: it passes fluid
//! without altering the topology of the network. Its hydraulic effect is captured
//! by a user-specified **insertion resistance** `R_insertion ≥ 0` [Pa·s/m³],
//! representing the pressure drop introduced by the sensor body (tubing connectors,
//! flow-cell housing, etc.):
//!
//! ```text
//! ΔP_sensor = R_insertion · Q
//! ```
//!
//! For a truly non-intrusive sensor `R_insertion = 0`; for a thermal mass-flow
//! sensor or Coriolis device, the insertion loss must be characterised
//! experimentally and provided by the user.
//!
//! ## Invariants
//!
//! - `resistance ≥ 0` (non-negative; zero is allowed for ideal sensors)
//! - `range > 0` (measurement range must be strictly positive)
//!
//! ## Reference
//!
//! Gravesen, P., Branebjerg, J., & Jensen, O. S. (1993). Microfluidics —
//! a review. *Journal of Micromechanics and Microengineering*, 3(4), 168–182.

use super::Component;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::Fluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Sensor type enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SensorType {
    /// Volumetric flow rate sensor
    Flow,
    /// Pressure sensor (zero insertion resistance)
    Pressure,
    /// Temperature sensor (zero insertion resistance)
    Temperature,
    /// Concentration / optical sensor
    Concentration,
}

/// Flow sensor component — ideal measurement device with user-specified insertion resistance.
///
/// See module documentation for the governing theorem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowSensor<T: RealField + Copy> {
    /// Insertion resistance [Pa·s/m³] (≥ 0; use 0.0 for ideal sensor)
    pub resistance: T,
    /// Measurement range [m³/s] (must be > 0)
    pub range: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> FlowSensor<T> {
    /// Create a new flow sensor with validated parameters.
    ///
    /// # Errors
    ///
    /// Returns `Error::InvalidConfiguration` if `resistance < 0` or `range ≤ 0`.
    pub fn new(resistance: T, range: T) -> Result<Self> {
        if resistance < T::zero() {
            return Err(Error::InvalidConfiguration(
                "FlowSensor resistance must be ≥ 0 (use 0.0 for ideal sensor)".into(),
            ));
        }
        if range <= T::zero() {
            return Err(Error::InvalidConfiguration(
                "FlowSensor measurement range must be > 0".into(),
            ));
        }
        Ok(Self {
            resistance,
            range,
            parameters: HashMap::new(),
        })
    }

    /// Return `true` if the measured flow rate exceeds the sensor range.
    #[must_use]
    pub fn is_overrange(&self, flow_rate: T) -> bool {
        flow_rate.abs() > self.range
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for FlowSensor<T> {
    /// Returns the sensor insertion resistance [Pa·s/m³].
    ///
    /// **Invariant**: `R_insertion ≥ 0`. A return value of `0.0` denotes
    /// an ideal non-intrusive sensor.
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        self.resistance
    }

    fn component_type(&self) -> &'static str {
        "FlowSensor"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "resistance" => {
                if value < T::zero() {
                    return Err(Error::InvalidConfiguration(
                        "FlowSensor resistance must be ≥ 0".into(),
                    ));
                }
                self.resistance = value;
            }
            "range" => {
                if value <= T::zero() {
                    return Err(Error::InvalidConfiguration(
                        "FlowSensor range must be > 0".into(),
                    ));
                }
                self.range = value;
            }
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }
}
