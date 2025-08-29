//! Component type definitions for network edges

use serde::{Deserialize, Serialize};

/// Physical component type for network edges
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum ComponentType {
    /// Standard pipe segment
    Pipe,
    /// Flow control valve
    Valve,
    /// Pump or compressor
    Pump,
    /// Junction or fitting
    Junction,
    /// Heat exchanger
    HeatExchanger,
    /// Flow meter or sensor
    Sensor,
    /// Reservoir or tank
    Reservoir,
    /// Generic component
    Generic,
}

impl ComponentType {
    /// Get string representation for display
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Pipe => "pipe",
            Self::Valve => "valve",
            Self::Pump => "pump",
            Self::Junction => "junction",
            Self::HeatExchanger => "heat_exchanger",
            Self::Sensor => "sensor",
            Self::Reservoir => "reservoir",
            Self::Generic => "generic",
        }
    }
}

impl Default for ComponentType {
    fn default() -> Self {
        Self::Pipe
    }
}
