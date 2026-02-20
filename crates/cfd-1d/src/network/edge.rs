//! Network edge definitions

use nalgebra::RealField;
use cfd_schematics::domain::model::EdgeKind;
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

use cfd_schematics::domain::model::ChannelSpec;
use num_traits::FromPrimitive;

impl<T: RealField + Copy + FromPrimitive> From<&ChannelSpec> for Edge<T> {
    fn from(spec: &ChannelSpec) -> Self {
        let mut resistance = T::from_f64(spec.resistance).unwrap_or(T::zero());
        let mut quad_coeff = T::from_f64(spec.quad_coeff).unwrap_or(T::zero());

        match spec.kind {
            EdgeKind::Valve => {
                if let Some(cv) = spec.valve_cv {
                    // Heuristic: k = 1/Cv^2 for initial validation
                    if cv > 0.0 {
                        let cv_t = T::from_f64(cv).unwrap_or(T::one());
                        quad_coeff = T::one() / (cv_t * cv_t);
                    }
                }
                // Ensure non-zero to pass builder validation if nothing else set
                if resistance.abs() < T::default_epsilon() && quad_coeff.abs() < T::default_epsilon() {
                     resistance = T::from_f64(1e-6).unwrap_or(T::zero());
                }
            }
            EdgeKind::Pump => {
                // Pumps need non-zero internal resistance for solver stability/builder validation
                if resistance.abs() < T::default_epsilon() && quad_coeff.abs() < T::default_epsilon() {
                     resistance = T::from_f64(1e-6).unwrap_or(T::zero());
                }
            }
            _ => {}
        }
        
        // Ensure resistance is non-zero for numerical stability (Newton-Raphson at Q=0)
        let min_r = T::from_f64(1e-8).unwrap_or(T::default_epsilon());
        if resistance.abs() < min_r {
            resistance = min_r;
        }

        Self {
            id: spec.id.as_str().to_string(),
            edge_type: spec.kind,
            flow_rate: T::zero(),
            resistance,
            quad_coeff,
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
