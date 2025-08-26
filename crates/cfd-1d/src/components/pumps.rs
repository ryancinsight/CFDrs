//! Pump components for microfluidic networks

use super::{constants, Component};
use cfd_core::{Fluid, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
/// Pump type enumeration
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum PumpType {
    /// Syringe pump
    Syringe,
    /// Peristaltic pump
    Peristaltic,
    /// Diaphragm pump
    Diaphragm,
    /// Electroosmotic pump
    Electroosmotic,
}
/// Micropump component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micropump<T: RealField + Copy> {
    /// Maximum flow rate [m³/s]
    pub max_flow_rate: T,
    /// Maximum pressure [Pa]
    pub max_pressure: T,
    /// Pump efficiency [-]
    pub efficiency: T,
    /// Current operating point (0-1)
    pub operating_point: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
impl<T: RealField + Copy + FromPrimitive + Float> Micropump<T> {
    /// Create a new micropump
    pub fn new(max_flow_rate: T, max_pressure: T) -> Self {
        Self {
            max_flow_rate,
            max_pressure,
            efficiency: T::from_f64(constants::DEFAULT_PUMP_EFFICIENCY).unwrap_or_else(T::one),
            operating_point: T::one(),
            parameters: HashMap::new(),
        }
    }
impl<T: RealField + Copy + FromPrimitive + Float> Component<T> for Micropump<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        // Pumps provide negative resistance (pressure source)
        -self.max_pressure / self.max_flow_rate
    fn component_type(&self) -> &str {
        "Micropump"
    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "max_flow_rate" => self.max_flow_rate = value,
            "max_pressure" => self.max_pressure = value,
            "efficiency" => self.efficiency = value,
            "operating_point" => self.operating_point = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        Ok(())
    fn is_active(&self) -> bool {
        true
