//! Microfluidic and millifluidic components for 1D network modeling.
//!
//! This module provides models for various microfluidic components including
//! channels, pumps, valves, mixers, and sensors with their characteristic
//! hydraulic properties and behaviors.

use cfd_core::error::Result;
use cfd_core::fluid::{ConstantPropertyFluid, Fluid};
use nalgebra::RealField;
use std::collections::HashMap;

// Re-export submodules
pub mod channels;
pub mod constants;
pub mod factory;
pub mod mixers;
pub mod pumps;
pub mod sensors;
pub mod valves;

// Re-export commonly used types
pub use channels::{CircularChannel, RectangularChannel};
pub use factory::ComponentFactory;
pub use mixers::{Micromixer, MixerType};
pub use pumps::{Micropump, PumpType};
pub use sensors::{FlowSensor, SensorType};
pub use valves::{Microvalve, ValveType};

/// Trait for all microfluidic components
pub trait Component<T: RealField + Copy> {
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
