//! Network edge definitions

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Edge types in the network
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EdgeType {
    /// Standard pipe
    Pipe,
    /// Control valve
    Valve,
    /// Pump for pressure boost
    Pump,
}

/// Edge in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge<T: RealField + Copy> {
    /// Unique identifier
    pub id: String,
    /// Type of edge
    pub edge_type: EdgeType,
    /// Flow rate through the edge
    pub flow_rate: T,
    /// Resistance coefficient
    pub resistance: T,
}

impl<T: RealField + Copy> Edge<T> {
    /// Create a new edge
    #[must_use]
    pub fn new(id: String, edge_type: EdgeType) -> Self {
        Self {
            id,
            edge_type,
            flow_rate: T::zero(),
            resistance: T::one(),
        }
    }
}

/// Channel-specific properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelProperties<T: RealField + Copy> {
    /// Length of the channel
    pub length: T,
    /// Diameter or characteristic dimension
    pub diameter: T,
    /// Surface roughness
    pub roughness: T,
}
