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
    /// Quadratic loss coefficient for dissipative term in the form ΔP = R·Q + k·Q|Q|
    /// where `k = quad_coeff ≥ 0`. This ensures the nonlinear loss opposes the flow
    /// direction and preserves monotonic positivity in the effective resistance
    /// linearization: R_eff = R + 2k|Q_k|.
    pub quad_coeff: T,
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
            quad_coeff: T::zero(),
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
