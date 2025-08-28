//! Network edge definitions

use crate::channel::ChannelGeometry;
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
    pub fn new(id: String, edge_type: EdgeType) -> Self {
        Self {
            id,
            edge_type,
            flow_rate: T::zero(),
            resistance: T::one(),
        }
    }
}

/// Properties for edges in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EdgeProperties<T: RealField + Copy> {
    /// Channel geometry
    pub geometry: Option<ChannelGeometry<T>>,
    /// Pressure drop across edge
    pub pressure_drop: T,
    /// Additional metadata
    pub metadata: std::collections::HashMap<String, T>,
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

impl<T: RealField + Copy> Default for EdgeProperties<T> {
    fn default() -> Self {
        Self {
            geometry: None,
            pressure_drop: T::zero(),
            metadata: std::collections::HashMap::new(),
        }
    }
}
