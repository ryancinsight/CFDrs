//! Network edge definitions

use cfd_schematics::domain::model::EdgeKind;
use nalgebra::RealField;
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
    /// Cross-sectional area \[m²\] derived from the channel geometry.
    /// Used for junction minor-loss K-factor corrections: ΔP = K·ρQ²/(2A²).
    pub area: T,
}

use cfd_schematics::domain::model::ChannelSpec;
use num_traits::FromPrimitive;

impl<T: RealField + Copy + FromPrimitive> From<&ChannelSpec> for Edge<T> {
    fn from(spec: &ChannelSpec) -> Self {
        let mut resistance =
            T::from_f64(spec.resistance).expect("Mathematical constant conversion compromised");
        let mut quad_coeff =
            T::from_f64(spec.quad_coeff).expect("Mathematical constant conversion compromised");

        match spec.kind {
            EdgeKind::Valve => {
                if let Some(cv) = spec.valve_cv {
                    // Use the Microvalve quadratic-law coefficient from the valve theorem.
                    if cv > 0.0 {
                        let cv_t =
                            T::from_f64(cv).expect("Mathematical constant conversion compromised");
                        quad_coeff = T::one() / (cv_t * cv_t);
                    }
                }
                // Ensure non-zero to pass builder validation if nothing else set
                if resistance.abs() < T::default_epsilon()
                    && quad_coeff.abs() < T::default_epsilon()
                {
                    resistance =
                        T::from_f64(1e-6).expect("Mathematical constant conversion compromised");
                }
            }
            EdgeKind::Pump | EdgeKind::Pipe => {}
        }

        let area = match spec.cross_section {
            cfd_schematics::domain::model::CrossSectionSpec::Circular { diameter_m } => {
                T::from_f64(std::f64::consts::PI * (diameter_m / 2.0).powi(2)).unwrap_or(T::zero())
            }
            cfd_schematics::domain::model::CrossSectionSpec::Rectangular { width_m, height_m } => {
                T::from_f64(width_m * height_m)
                    .expect("Mathematical constant conversion compromised")
            }
        };

        Self {
            id: spec.id.as_str().to_string(),
            edge_type: spec.kind,
            flow_rate: T::zero(),
            resistance,
            quad_coeff,
            area,
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
            area: T::zero(),
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
