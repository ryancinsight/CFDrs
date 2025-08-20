//! Microfluidic and millifluidic components for 1D network modeling.
//!
//! This module provides models for various microfluidic components including
//! channels, pumps, valves, mixers, and sensors with their characteristic
//! hydraulic properties and behaviors.

use cfd_core::{Error, Result, Fluid};
use nalgebra::RealField;
use std::collections::HashMap;

// Re-export submodules
pub mod channels;
pub mod pumps;
pub mod valves;
pub mod sensors;
pub mod mixers;
pub mod factory;

// Re-export commonly used types
pub use channels::{RectangularChannel, CircularChannel};
pub use pumps::Micropump;
pub use valves::Microvalve;
pub use sensors::FlowSensor;
pub use mixers::Micromixer;
pub use factory::ComponentFactory;

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

/// Constants for component physics
pub mod constants {
    /// Default surface roughness for smooth channels [m]
    pub const DEFAULT_ROUGHNESS: f64 = 1e-6;
    
    /// Minimum Reynolds number for laminar flow
    pub const RE_LAMINAR_MIN: f64 = 0.1;
    
    /// Maximum Reynolds number for laminar flow
    pub const RE_LAMINAR_MAX: f64 = 2300.0;
    
    /// Transition Reynolds number
    pub const RE_TRANSITION: f64 = 2300.0;
    
    /// Maximum Reynolds number for transition flow
    pub const RE_TURBULENT_MIN: f64 = 4000.0;
    
    /// Default pump efficiency
    pub const DEFAULT_PUMP_EFFICIENCY: f64 = 0.7;
    
    /// Default valve flow coefficient
    pub const DEFAULT_VALVE_CV: f64 = 0.1;
    
    /// Default mixing efficiency
    pub const DEFAULT_MIXING_EFFICIENCY: f64 = 0.95;
}