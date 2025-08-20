//! Microfluidic and millifluidic components for 1D network modeling.
//!
//! This module provides models for various microfluidic components including
//! channels, pumps, valves, mixers, and sensors with their characteristic
//! hydraulic properties and behaviors.

use cfd_core::{Error, Result, Fluid};
use nalgebra::{RealField, ComplexField};
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Trait for all microfluidic components
pub trait Component<T: RealField> {
    /// Get the hydraulic resistance of the component
    fn resistance(&self, fluid: &Fluid<T>) -> T;

    /// Get the pressure drop across the component for a given flow rate
    fn pressure_drop(&self, flow_rate: T, fluid: &Fluid<T>) -> T {
        flow_rate * self.resistance(fluid)
    }

    /// Get the component type identifier
    fn component_type(&self) -> &str;

    /// Get component parameters
    fn parameters(&self) -> &HashMap<String, T>;

    /// Update component parameters
    fn set_parameter(&mut self, key: &str, value: T) -> Result<()>;

    /// Check if the component is active (e.g., pump, valve)
    fn is_active(&self) -> bool {
        false
    }

    /// Get the volume of the component (for transient analysis)
    fn volume(&self) -> Option<T> {
        None
    }
}

/// Rectangular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RectangularChannel<T: RealField> {
    /// Channel length [m]
    pub length: T,
    /// Channel width [m]
    pub width: T,
    /// Channel height [m]
    pub height: T,
    /// Surface roughness [m]
    pub roughness: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> RectangularChannel<T> {
    /// Create a new rectangular channel
    pub fn new(length: T, width: T, height: T, roughness: T) -> Self {
        Self {
            length,
            width,
            height,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Create a square channel
    pub fn square(length: T, side: T, roughness: T) -> Self {
        Self::new(length, side.clone(), side, roughness)
    }

    /// Get cross-sectional area
    pub fn area(&self) -> T {
        self.width.clone() * self.height.clone()
    }

    /// Get hydraulic diameter
    pub fn hydraulic_diameter(&self) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        two * self.area() / (self.width.clone() + self.height.clone())
    }

    /// Get aspect ratio (width/height)
    pub fn aspect_ratio(&self) -> T {
        self.width.clone() / self.height.clone()
    }

    /// Calculate friction factor for laminar flow
    fn friction_factor_laminar(&self) -> T {
        let alpha = self.aspect_ratio();
        let one = T::one();
        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        let twentyfour = T::from_f64(24.0).unwrap_or_else(|| T::zero());
        let sixtyfour = T::from_f64(64.0).unwrap_or_else(|| T::zero());

        if alpha >= one {
            // Wide channel approximation
            twentyfour * (one.clone() - T::from_f64(1.3553).unwrap_or_else(|| T::zero()) * alpha.clone() +
                         T::from_f64(1.9467).unwrap_or_else(|| T::zero()) * alpha.clone() * alpha.clone() -
                         T::from_f64(1.7012).unwrap_or_else(|| T::zero()) * ComplexField::powf(alpha.clone(), three.clone()) +
                         T::from_f64(0.9564).unwrap_or_else(|| T::zero()) * ComplexField::powf(alpha.clone(), four.clone()) -
                         T::from_f64(0.2537).unwrap_or_else(|| T::zero()) * ComplexField::powf(alpha, T::from_f64(5.0).unwrap_or_else(|| T::zero())))
        } else {
            // Tall channel (alpha < 1)
            let inv_alpha = one / alpha;
            sixtyfour / inv_alpha.clone() - T::from_f64(13.74).unwrap_or_else(|| T::zero()) +
            T::from_f64(7.37).unwrap_or_else(|| T::zero()) * inv_alpha.clone() -
            T::from_f64(1.27).unwrap_or_else(|| T::zero()) * inv_alpha.clone() * inv_alpha
        }
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for RectangularChannel<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        let dh = self.hydraulic_diameter();
        let area = self.area();
        let f = self.friction_factor_laminar();

        // Hydraulic resistance: R = (f * L * μ) / (ρ * A * Dh^2)
        // For pressure-driven flow: R = (f * L) / (A * Dh^2) * (μ/ρ) = (f * L * ν) / (A * Dh^2)
        // Simplified for laminar flow: R = (128 * μ * L) / (π * D^4) for circular
        // For rectangular: R = (f * L * μ) / (A * Dh^2 * ρ)

        // Use actual operating temperature instead of T::zero()
        let temperature = T::from_f64(293.15).unwrap_or_else(|| T::zero()); // Default to 20°C if not specified
        let kinematic_viscosity = fluid.dynamic_viscosity(temperature)
            .unwrap_or_else(|_| T::from_f64(0.001).unwrap_or_else(|| T::one())) / fluid.density;
        let resistance = f * self.length.clone() * kinematic_viscosity / (area * dh.clone() * dh);

        // Ensure positive resistance
        if resistance > T::zero() {
            resistance
        } else {
            T::from_f64(1e-12).unwrap_or_else(|| T::zero()) // Minimum resistance
        }
    }

    fn component_type(&self) -> &str {
        "RectangularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "width" => self.width = value,
            "height" => self.height = value,
            "roughness" => self.roughness = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.length.clone() * self.area())
    }
}

/// Circular microchannel component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CircularChannel<T: RealField> {
    /// Channel length [m]
    pub length: T,
    /// Channel diameter [m]
    pub diameter: T,
    /// Surface roughness [m]
    pub roughness: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + FromPrimitive + num_traits::Float> CircularChannel<T> {
    /// Create a new circular channel
    pub fn new(length: T, diameter: T, roughness: T) -> Self {
        Self {
            length,
            diameter,
            roughness,
            parameters: HashMap::new(),
        }
    }

    /// Get cross-sectional area
    pub fn area(&self) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        pi * self.diameter.clone() * self.diameter.clone() / four
    }

    /// Get radius
    pub fn radius(&self) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        self.diameter.clone() / two
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for CircularChannel<T> {
    fn resistance(&self, fluid: &Fluid<T>) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero());
        let onehundredtwentyeight = T::from_f64(128.0).unwrap_or_else(|| T::zero());

        // Hagen-Poiseuille equation: R = (128 * μ * L) / (π * D^4)
        let d4 = ComplexField::powf(self.diameter.clone(), T::from_f64(4.0).unwrap_or_else(|| T::zero()));
        // Use actual operating temperature instead of T::zero()
        let temperature = T::from_f64(293.15).unwrap_or_else(|| T::zero()); // Default to 20°C if not specified
        onehundredtwentyeight * fluid.dynamic_viscosity(temperature) * self.length.clone() / (pi * d4)
    }

    fn component_type(&self) -> &str {
        "CircularChannel"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "length" => self.length = value,
            "diameter" => self.diameter = value,
            "roughness" => self.roughness = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.area() * self.length.clone())
    }
}

/// Micropump component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micropump<T: RealField> {
    /// Maximum pressure rise [Pa]
    pub max_pressure: T,
    /// Maximum flow rate [m³/s]
    pub max_flow_rate: T,
    /// Pump efficiency [0-1]
    pub efficiency: T,
    /// Current operating point (0-1, where 1 = max)
    pub operating_point: T,
    /// Pump type
    pub pump_type: PumpType,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

/// Types of micropumps
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum PumpType {
    /// Peristaltic pump
    Peristaltic,
    /// Syringe pump
    Syringe,
    /// Electroosmotic pump
    Electroosmotic,
    /// Pressure-driven pump
    PressureDriven,
    /// Centrifugal pump
    Centrifugal,
}

impl<T: RealField + FromPrimitive + num_traits::Float> Micropump<T> {
    /// Create a new micropump
    pub fn new(max_pressure: T, max_flow_rate: T, pump_type: PumpType) -> Self {
        Self {
            max_pressure,
            max_flow_rate,
            efficiency: T::from_f64(0.8).unwrap_or_else(|| T::zero()),
            operating_point: T::from_f64(0.5).unwrap_or_else(|| T::zero()),
            pump_type,
            parameters: HashMap::new(),
        }
    }

    /// Create a syringe pump
    pub fn syringe_pump(max_pressure: T, max_flow_rate: T) -> Self {
        Self::new(max_pressure, max_flow_rate, PumpType::Syringe)
    }

    /// Create a peristaltic pump
    pub fn peristaltic_pump(max_pressure: T, max_flow_rate: T) -> Self {
        Self::new(max_pressure, max_flow_rate, PumpType::Peristaltic)
    }

    /// Get current pressure rise based on operating point
    pub fn current_pressure_rise(&self) -> T {
        self.max_pressure.clone() * self.operating_point.clone() * self.efficiency.clone()
    }

    /// Get current flow rate capacity
    pub fn current_flow_capacity(&self) -> T {
        self.max_flow_rate.clone() * self.operating_point.clone()
    }

    /// Set operating point (0-1)
    pub fn set_operating_point(&mut self, point: T) -> Result<()> {
        if point < T::zero() || point > T::one() {
            return Err(Error::InvalidConfiguration(
                "Operating point must be between 0 and 1".to_string()
            ));
        }
        self.operating_point = point;
        Ok(())
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for Micropump<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        // Pumps typically have very low resistance
        T::from_f64(1e-6).unwrap_or_else(|| T::zero())
    }

    fn pressure_drop(&self, _flow_rate: T, _fluid: &Fluid<T>) -> T {
        // Pumps provide pressure rise, not drop
        -self.current_pressure_rise()
    }

    fn component_type(&self) -> &str {
        "Micropump"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "max_pressure" => self.max_pressure = value,
            "max_flow_rate" => self.max_flow_rate = value,
            "efficiency" => {
                if value < T::zero() || value > T::one() {
                    return Err(Error::InvalidConfiguration(
                        "Efficiency must be between 0 and 1".to_string()
                    ));
                }
                self.efficiency = value;
            },
            "operating_point" => self.set_operating_point(value)?,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn is_active(&self) -> bool {
        true
    }
}

/// Microvalve component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Microvalve<T: RealField> {
    /// Base resistance when fully open [Pa·s/m³]
    pub base_resistance: T,
    /// Current opening fraction [0-1]
    pub opening: T,
    /// Valve type
    pub valve_type: ValveType,
    /// Response time [s]
    pub response_time: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

/// Types of microvalves
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum ValveType {
    /// Normally open valve
    NormallyOpen,
    /// Normally closed valve
    NormallyClosed,
    /// Proportional valve
    Proportional,
    /// On/off valve
    OnOff,
    /// Check valve (one-way)
    CheckValve,
}

impl<T: RealField + FromPrimitive + num_traits::Float> Microvalve<T> {
    /// Create a new microvalve
    pub fn new(base_resistance: T, valve_type: ValveType) -> Self {
        let initial_opening = match valve_type {
            ValveType::NormallyOpen => T::one(),
            ValveType::NormallyClosed => T::zero(),
            _ => T::from_f64(0.5).unwrap_or_else(|| T::zero()),
        };

        Self {
            base_resistance,
            opening: initial_opening,
            valve_type,
            response_time: T::from_f64(0.001).unwrap_or_else(|| T::zero()), // 1 ms default
            parameters: HashMap::new(),
        }
    }

    /// Create a normally open valve
    pub fn normally_open(base_resistance: T) -> Self {
        Self::new(base_resistance, ValveType::NormallyOpen)
    }

    /// Create a normally closed valve
    pub fn normally_closed(base_resistance: T) -> Self {
        Self::new(base_resistance, ValveType::NormallyClosed)
    }

    /// Set valve opening (0 = closed, 1 = fully open)
    pub fn set_opening(&mut self, opening: T) -> Result<()> {
        if opening < T::zero() || opening > T::one() {
            return Err(Error::InvalidConfiguration(
                "Valve opening must be between 0 and 1".to_string()
            ));
        }
        self.opening = opening;
        Ok(())
    }

    /// Open the valve
    pub fn open(&mut self) {
        self.opening = T::one();
    }

    /// Close the valve
    pub fn close(&mut self) {
        self.opening = T::zero();
    }

    /// Check if valve is open
    pub fn is_open(&self) -> bool {
        self.opening > T::zero()
    }

    /// Check if valve is closed
    pub fn is_closed(&self) -> bool {
        self.opening == T::zero()
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for Microvalve<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        if self.opening == T::zero() {
            // Closed valve - very high resistance
            T::from_f64(1e12).unwrap_or_else(|| T::zero())
        } else {
            // Resistance inversely proportional to opening
            self.base_resistance.clone() / self.opening.clone()
        }
    }

    fn component_type(&self) -> &str {
        "Microvalve"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "base_resistance" => self.base_resistance = value,
            "opening" => self.set_opening(value)?,
            "response_time" => self.response_time = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn is_active(&self) -> bool {
        true
    }
}

/// Flow sensor component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlowSensor<T: RealField> {
    /// Sensor accuracy [0-1]
    pub accuracy: T,
    /// Measurement range [m³/s]
    pub range: T,
    /// Response time [s]
    pub response_time: T,
    /// Sensor type
    pub sensor_type: SensorType,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

/// Types of flow sensors
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SensorType {
    /// Thermal flow sensor
    Thermal,
    /// Pressure differential sensor
    PressureDifferential,
    /// Ultrasonic sensor
    Ultrasonic,
    /// Optical sensor
    Optical,
    /// Coriolis sensor
    Coriolis,
}

impl<T: RealField + FromPrimitive + num_traits::Float> FlowSensor<T> {
    /// Create a new flow sensor
    pub fn new(range: T, accuracy: T, sensor_type: SensorType) -> Self {
        Self {
            accuracy,
            range,
            response_time: T::from_f64(0.001).unwrap_or_else(|| T::zero()), // 1 ms default
            sensor_type,
            parameters: HashMap::new(),
        }
    }

    /// Create a thermal flow sensor
    pub fn thermal(range: T, accuracy: T) -> Self {
        Self::new(range, accuracy, SensorType::Thermal)
    }

    /// Simulate sensor measurement with noise
    pub fn measure(&self, true_flow_rate: T) -> T {
        // Add measurement noise based on accuracy
        let noise_factor = T::one() - self.accuracy.clone();
        let noise = (T::from_f64(rand::random::<f64>()).expect("FIXME: Add proper error handling") - T::from_f64(0.5).unwrap_or_else(|| T::zero())) * noise_factor;
        true_flow_rate * (T::one() + noise)
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for FlowSensor<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        // Sensors should have minimal resistance
        T::from_f64(1e-9).unwrap_or_else(|| T::zero())
    }

    fn component_type(&self) -> &str {
        "FlowSensor"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "accuracy" => {
                if value < T::zero() || value > T::one() {
                    return Err(Error::InvalidConfiguration(
                        "Accuracy must be between 0 and 1".to_string()
                    ));
                }
                self.accuracy = value;
            },
            "range" => self.range = value,
            "response_time" => self.response_time = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }
}

/// Micromixer component
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micromixer<T: RealField> {
    /// Mixing chamber volume [m³]
    pub volume: T,
    /// Mixing efficiency [0-1]
    pub efficiency: T,
    /// Pressure drop coefficient
    pub pressure_drop_coeff: T,
    /// Mixer type
    pub mixer_type: MixerType,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

/// Types of micromixers
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum MixerType {
    /// T-junction mixer
    TJunction,
    /// Y-junction mixer
    YJunction,
    /// Serpentine mixer
    Serpentine,
    /// Chaotic mixer
    Chaotic,
    /// Active mixer (with external energy)
    Active,
}

impl<T: RealField + FromPrimitive + num_traits::Float> Micromixer<T> {
    /// Create a new micromixer
    pub fn new(volume: T, efficiency: T, mixer_type: MixerType) -> Self {
        Self {
            volume,
            efficiency,
            pressure_drop_coeff: T::from_f64(1.0).unwrap_or_else(|| T::zero()),
            mixer_type,
            parameters: HashMap::new(),
        }
    }

    /// Create a T-junction mixer
    pub fn t_junction(volume: T, efficiency: T) -> Self {
        Self::new(volume, efficiency, MixerType::TJunction)
    }

    /// Create a serpentine mixer
    pub fn serpentine(volume: T, efficiency: T) -> Self {
        let mut mixer = Self::new(volume, efficiency, MixerType::Serpentine);
        mixer.pressure_drop_coeff = T::from_f64(2.0).unwrap_or_else(|| T::zero()); // Higher pressure drop
        mixer
    }

    /// Calculate mixing time based on flow rates
    pub fn mixing_time(&self, total_flow_rate: T) -> T {
        if total_flow_rate == T::zero() {
            T::from_f64(f64::INFINITY).unwrap_or_else(|| T::zero())
        } else {
            self.volume.clone() / total_flow_rate
        }
    }
}

impl<T: RealField + FromPrimitive + num_traits::Float> Component<T> for Micromixer<T> {
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        // Resistance based on volume and pressure drop coefficient
        let base_resistance = T::from_f64(1e6).unwrap_or_else(|| T::zero()); // Base value
        base_resistance * self.pressure_drop_coeff.clone() / self.volume.clone()
    }

    fn component_type(&self) -> &str {
        "Micromixer"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "volume" => self.volume = value,
            "efficiency" => {
                if value < T::zero() || value > T::one() {
                    return Err(Error::InvalidConfiguration(
                        "Efficiency must be between 0 and 1".to_string()
                    ));
                }
                self.efficiency = value;
            },
            "pressure_drop_coeff" => self.pressure_drop_coeff = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn volume(&self) -> Option<T> {
        Some(self.volume.clone())
    }
}

/// Component factory for creating standard components
pub struct ComponentFactory;

impl ComponentFactory {
    /// Create a standard rectangular microchannel
    pub fn standard_microchannel<T: RealField + FromPrimitive + num_traits::Float>(
        length: T,
        width: T,
        height: T,
    ) -> RectangularChannel<T> {
        RectangularChannel::new(length, width, height, T::from_f64(1e-6).unwrap_or_else(|| T::zero()))
    }

    /// Create a standard syringe pump
    pub fn standard_syringe_pump<T: RealField + FromPrimitive + num_traits::Float>() -> Micropump<T> {
        Micropump::syringe_pump(
            T::from_f64(1e5).unwrap_or_else(|| T::zero()),  // 100 kPa max pressure
            T::from_f64(1e-6).unwrap_or_else(|| T::zero()), // 1 μL/s max flow rate
        )
    }

    /// Create a standard microvalve
    pub fn standard_microvalve<T: RealField + FromPrimitive + num_traits::Float>() -> Microvalve<T> {
        Microvalve::normally_closed(T::from_f64(1e8).unwrap_or_else(|| T::zero())) // 100 MPa·s/m³
    }

    /// Create a standard flow sensor
    pub fn standard_flow_sensor<T: RealField + FromPrimitive + num_traits::Float>() -> FlowSensor<T> {
        FlowSensor::thermal(
            T::from_f64(1e-5).unwrap_or_else(|| T::zero()), // 10 μL/s range
            T::from_f64(0.95).unwrap_or_else(|| T::zero()), // 95% accuracy
        )
    }

    /// Create a standard T-junction mixer
    pub fn standard_t_mixer<T: RealField + FromPrimitive + num_traits::Float>() -> Micromixer<T> {
        Micromixer::t_junction(
            T::from_f64(1e-9).unwrap_or_else(|| T::zero()), // 1 nL volume
            T::from_f64(0.8).unwrap_or_else(|| T::zero()),  // 80% efficiency
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rectangular_channel() {
        let channel = RectangularChannel::new(0.001, 100e-6, 50e-6, 1e-6);

        assert_relative_eq!(channel.area(), 5e-9, epsilon = 1e-15);
        assert_relative_eq!(channel.hydraulic_diameter(), 200e-6 / 3.0, epsilon = 1e-15);
        assert_relative_eq!(channel.aspect_ratio(), 2.0, epsilon = 1e-10);

        let fluid = Fluid::water();
        let resistance = channel.resistance(&fluid);
        assert!(resistance > 0.0);

        assert_eq!(channel.component_type(), "RectangularChannel");
        assert!(channel.volume().is_some());
    }

    #[test]
    fn test_square_channel() {
        let channel = RectangularChannel::square(0.001, 100e-6, 1e-6);

        assert_relative_eq!(channel.aspect_ratio(), 1.0, epsilon = 1e-10);
        assert_relative_eq!(channel.area(), 1e-8, epsilon = 1e-15);
    }

    #[test]
    fn test_circular_channel() {
        let channel = CircularChannel::new(0.001, 100e-6, 1e-6);

        let expected_area = std::f64::consts::PI * (50e-6_f64).powi(2);
        assert_relative_eq!(channel.area(), expected_area, epsilon = 1e-15);
        assert_relative_eq!(channel.radius(), 50e-6, epsilon = 1e-15);

        let fluid = Fluid::water();
        let resistance = channel.resistance(&fluid);
        assert!(resistance > 0.0);

        assert_eq!(channel.component_type(), "CircularChannel");
        assert!(channel.volume().is_some());
    }

    #[test]
    fn test_micropump() {
        let mut pump = Micropump::syringe_pump(1000.0, 1e-6);

        assert_eq!(pump.pump_type, PumpType::Syringe);
        assert!(pump.is_active());

        // Test operating point
        pump.set_operating_point(0.8).expect("FIXME: Add proper error handling");
        assert_relative_eq!(pump.operating_point, 0.8, epsilon = 1e-10);

        // Test pressure rise calculation
        let pressure_rise = pump.current_pressure_rise();
        assert_relative_eq!(pressure_rise, 1000.0 * 0.8 * 0.8, epsilon = 1e-10); // max * operating * efficiency

        // Test invalid operating point
        assert!(pump.set_operating_point(1.5).is_err());
        assert!(pump.set_operating_point(-0.1).is_err());

        let fluid = Fluid::water();
        assert!(pump.resistance(&fluid) > 0.0);
        assert_eq!(pump.component_type(), "Micropump");
    }

    #[test]
    fn test_microvalve() {
        let mut valve = Microvalve::normally_closed(1e8);

        assert_eq!(valve.valve_type, ValveType::NormallyClosed);
        assert!(valve.is_closed());
        assert!(!valve.is_open());
        assert!(valve.is_active());

        // Test opening the valve
        valve.open();
        assert!(valve.is_open());
        assert!(!valve.is_closed());

        // Test partial opening
        valve.set_opening(0.5).expect("FIXME: Add proper error handling");
        assert_relative_eq!(valve.opening, 0.5, epsilon = 1e-10);

        // Test resistance calculation
        let fluid = Fluid::water();
        let resistance_half_open = valve.resistance(&fluid);
        assert_relative_eq!(resistance_half_open, 2e8, epsilon = 1e-6); // base_resistance / opening

        // Test closed valve resistance
        valve.close();
        let resistance_closed = valve.resistance(&fluid);
        assert!(resistance_closed > 1e10); // Very high resistance when closed

        // Test invalid opening
        assert!(valve.set_opening(1.5).is_err());
        assert!(valve.set_opening(-0.1).is_err());

        assert_eq!(valve.component_type(), "Microvalve");
    }

    #[test]
    fn test_normally_open_valve() {
        let valve = Microvalve::normally_open(1e6);

        assert_eq!(valve.valve_type, ValveType::NormallyOpen);
        assert!(valve.is_open());
        assert!(!valve.is_closed());
        assert_relative_eq!(valve.opening, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_flow_sensor() {
        let mut sensor = FlowSensor::thermal(1e-5, 0.95);

        assert_eq!(sensor.sensor_type, SensorType::Thermal);
        assert_relative_eq!(sensor.accuracy, 0.95, epsilon = 1e-10);
        assert_relative_eq!(sensor.range, 1e-5, epsilon = 1e-15);

        // Test measurement (with some tolerance for noise)
        let true_flow = 5e-6;
        let measured_flow = sensor.measure(true_flow);
        let error: f64 = (measured_flow - true_flow).abs() / true_flow;
        assert!(error < 0.1); // Should be within 10% due to 95% accuracy

        // Test parameter setting
        sensor.set_parameter("accuracy", 0.99).expect("FIXME: Add proper error handling");
        assert_relative_eq!(sensor.accuracy, 0.99, epsilon = 1e-10);

        // Test invalid accuracy
        assert!(sensor.set_parameter("accuracy", 1.5).is_err());
        assert!(sensor.set_parameter("accuracy", -0.1).is_err());

        let fluid = Fluid::water();
        let resistance = sensor.resistance(&fluid);
        assert!(resistance < 1e-6); // Very low resistance

        assert_eq!(sensor.component_type(), "FlowSensor");
    }

    #[test]
    fn test_micromixer() {
        let mut mixer = Micromixer::t_junction(1e-9, 0.8);

        assert_eq!(mixer.mixer_type, MixerType::TJunction);
        assert_relative_eq!(mixer.efficiency, 0.8, epsilon = 1e-10);
        assert_relative_eq!(mixer.volume, 1e-9, epsilon = 1e-15);

        // Test mixing time calculation
        let flow_rate = 1e-6;
        let mixing_time = mixer.mixing_time(flow_rate);
        assert_relative_eq!(mixing_time, 1e-3, epsilon = 1e-10); // volume / flow_rate

        // Test zero flow rate
        let mixing_time_zero: f64 = mixer.mixing_time(0.0);
        assert!(mixing_time_zero.is_infinite());

        // Test parameter setting
        mixer.set_parameter("efficiency", 0.9).expect("FIXME: Add proper error handling");
        assert_relative_eq!(mixer.efficiency, 0.9, epsilon = 1e-10);

        // Test invalid efficiency
        assert!(mixer.set_parameter("efficiency", 1.5).is_err());
        assert!(mixer.set_parameter("efficiency", -0.1).is_err());

        let fluid = Fluid::water();
        let resistance = mixer.resistance(&fluid);
        assert!(resistance > 0.0);

        assert_eq!(mixer.component_type(), "Micromixer");
        assert!(mixer.volume().is_some());
    }

    #[test]
    fn test_serpentine_mixer() {
        let mixer = Micromixer::serpentine(2e-9, 0.9);

        assert_eq!(mixer.mixer_type, MixerType::Serpentine);
        assert_relative_eq!(mixer.pressure_drop_coeff, 2.0, epsilon = 1e-10); // Higher than T-junction
    }

    #[test]
    fn test_component_factory() {
        let channel = ComponentFactory::standard_microchannel(0.001, 100e-6, 50e-6);
        assert_eq!(channel.component_type(), "RectangularChannel");

        let pump = ComponentFactory::standard_syringe_pump::<f64>();
        assert_eq!(pump.component_type(), "Micropump");
        assert_eq!(pump.pump_type, PumpType::Syringe);

        let valve = ComponentFactory::standard_microvalve::<f64>();
        assert_eq!(valve.component_type(), "Microvalve");
        assert_eq!(valve.valve_type, ValveType::NormallyClosed);

        let sensor = ComponentFactory::standard_flow_sensor::<f64>();
        assert_eq!(sensor.component_type(), "FlowSensor");
        assert_eq!(sensor.sensor_type, SensorType::Thermal);

        let mixer = ComponentFactory::standard_t_mixer::<f64>();
        assert_eq!(mixer.component_type(), "Micromixer");
        assert_eq!(mixer.mixer_type, MixerType::TJunction);
    }

    #[test]
    fn test_component_parameters() {
        let mut channel = RectangularChannel::new(0.001, 100e-6, 50e-6, 1e-6);

        // Test setting built-in parameters
        channel.set_parameter("length", 0.002).expect("FIXME: Add proper error handling");
        assert_relative_eq!(channel.length, 0.002, epsilon = 1e-15);

        // Test setting custom parameters
        channel.set_parameter("custom_param", 42.0).expect("FIXME: Add proper error handling");
        assert_eq!(channel.parameters().get("custom_param"), Some(&42.0));

        // Test getting parameters
        let params = channel.parameters();
        assert!(params.contains_key("custom_param"));
    }

    #[test]
    fn test_pressure_drop_calculation() {
        let channel = CircularChannel::new(0.001, 100e-6, 1e-6);
        let fluid = Fluid::water();

        let flow_rate = 1e-6; // 1 μL/s
        let pressure_drop = channel.pressure_drop(flow_rate, &fluid);
        let expected = flow_rate * channel.resistance(&fluid);

        assert_relative_eq!(pressure_drop, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_pump_pressure_drop() {
        let pump = Micropump::syringe_pump(1000.0, 1e-6);
        let fluid = Fluid::water();

        // Pump should provide negative pressure drop (pressure rise)
        let pressure_drop: f64 = pump.pressure_drop(1e-6, &fluid);
        assert!(pressure_drop < 0.0);
        assert_relative_eq!(pressure_drop.abs(), pump.current_pressure_rise(), epsilon = 1e-10);
    }
}
