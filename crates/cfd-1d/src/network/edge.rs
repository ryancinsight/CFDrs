//! Network edge definitions

use nalgebra::RealField;
use scheme::domain::model::EdgeKind;
use serde::{Deserialize, Serialize};

/// Edge in the network
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Edge<T: RealField + Copy> {
    /// Unique identifier
    pub id: String,
    /// Type of edge
    pub edge_type: EdgeKind,
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

use scheme::domain::model::ChannelSpec;
use num_traits::FromPrimitive;

impl<T: RealField + Copy + FromPrimitive> From<&ChannelSpec> for Edge<T> {
    fn from(spec: &ChannelSpec) -> Self {
        Self {
            id: spec.id.as_str().to_string(),
            edge_type: spec.kind,
            flow_rate: T::zero(),
            resistance: T::from_f64(spec.resistance).unwrap_or(T::one()),
            quad_coeff: T::from_f64(spec.quad_coeff).unwrap_or(T::zero()),
        }
    }
}

impl<T: RealField + Copy> Edge<T> {
    /// Create a new edge
    #[must_use]
    pub fn new(id: String, edge_type: EdgeKind) -> Self {
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
