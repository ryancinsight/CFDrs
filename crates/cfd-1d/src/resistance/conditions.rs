//! Flow conditions for resistance calculations

use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
/// Flow conditions for resistance calculations
#[derive(Debug, Clone)]
pub struct FlowConditions<T: RealField + Copy> {
    /// Reynolds number
    pub reynolds_number: Option<T>,
    /// Flow velocity [m/s]
    pub velocity: Option<T>,
    /// Flow rate [m³/s]
    pub flow_rate: Option<T>,
    /// Temperature [K]
    pub temperature: T,
    /// Pressure [Pa]
    pub pressure: T,
}
impl<T: RealField + Copy + FromPrimitive> FlowConditions<T> {
    /// Create default flow conditions
    pub fn default() -> Self {
        use cfd_core::constants::physical::{pressure, temperature};
        Self {
            reynolds_number: None,
            velocity: None,
            flow_rate: None,
            temperature: T::from_f64(temperature::CELSIUS_TO_KELVIN_OFFSET + 20.0)
                .unwrap_or_else(|| T::one()), // 20°C
            pressure: T::from_f64(pressure::STANDARD_ATMOSPHERE_PA).unwrap_or_else(|| T::one()), // 1 atm
        }
    }
    /// Set Reynolds number
    pub fn with_reynolds(mut self, re: T) -> Self {
        self.reynolds_number = Some(re);
        self
    /// Set flow velocity
    pub fn with_velocity(mut self, velocity: T) -> Self {
        self.velocity = Some(velocity);
    /// Set flow rate
    pub fn with_flow_rate(mut self, flow_rate: T) -> Self {
        self.flow_rate = Some(flow_rate);
