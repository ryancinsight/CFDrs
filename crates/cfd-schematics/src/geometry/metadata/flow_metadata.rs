use aequitas::systems::si::quantities::{Length, Pressure, VolumetricFlowRate};
use serde::{Deserialize, Serialize};

/// Fluid-flow metadata computed by the reduced-order solver.
#[derive(Debug, Clone, PartialEq)]
pub struct FlowMetadata {
    /// Volumetric flow rate.
    pub flow_rate: f64,
    /// Pressure drop across the channel.
    pub pressure_drop: f64,
    /// Dimensionless Reynolds number of the flow.
    pub reynolds_number: f64,
    /// Mean flow velocity.
    pub velocity: f64,
}

crate::impl_metadata!(FlowMetadata, "FlowMetadata");

/// Thermal metadata associated with a channel.
#[derive(Debug, Clone, PartialEq)]
pub struct ThermalMetadata {
    /// Channel temperature.
    pub temperature: f64,
    /// Convective heat transfer coefficient.
    pub heat_transfer_coefficient: f64,
    /// Material thermal conductivity.
    pub thermal_conductivity: f64,
}

crate::impl_metadata!(ThermalMetadata, "ThermalMetadata");

/// Manufacturing metadata associated with a channel.
#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturingMetadata {
    /// Toleranced channel width.
    pub width_tolerance: f64,
    /// Toleranced channel height.
    pub height_tolerance: f64,
    /// Surface roughness value.
    pub surface_roughness: f64,
    /// Manufacturing method identifier.
    pub manufacturing_method: String,
}

crate::impl_metadata!(ManufacturingMetadata, "ManufacturingMetadata");

/// Channel geometry metadata used by the reduced-order solver.
#[derive(Debug, Clone, PartialEq)]
pub struct ChannelGeometryMetadata {
    /// Equivalent channel diameter.
    pub channel_diameter_m: Length<f64>,
}

crate::impl_metadata!(ChannelGeometryMetadata, "ChannelGeometryMetadata");

/// Branch-level boundary condition metadata authored in schematics.
///
/// This metadata lets a blueprint mark an endpoint branch with an explicit
/// pressure or flow boundary condition instead of relying only on the legacy
/// `NodeKind` inference in reduced-order solvers.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BranchBoundarySpecification {
    /// Fixed pressure boundary.
    Pressure {
        /// Applied pressure at the branch.
        pressure_pa: Pressure<f64>,
    },
    /// Fixed volumetric flow boundary.
    ///
    /// The sign convention matches the 1D solver: positive values act as a
    /// source term, negative values as a sink term.
    FlowRate {
        /// Applied volumetric flow rate.
        flow_rate_m3_s: VolumetricFlowRate<f64>,
    },
}

/// Node-level branch boundary metadata.
///
/// Attach this to a `NodeSpec` in the schematics blueprint when a branch end
/// must override the default inlet/outlet pressure inference.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BranchBoundaryMetadata {
    /// Boundary condition specification applied at this branch end.
    pub boundary: BranchBoundarySpecification,
}

impl BranchBoundaryMetadata {
    /// Create a fixed-pressure branch boundary.
    #[must_use]
    pub fn pressure(pressure_pa: Pressure<f64>) -> Self {
        Self {
            boundary: BranchBoundarySpecification::Pressure { pressure_pa },
        }
    }

    /// Create a fixed-flow branch boundary.
    #[must_use]
    pub fn flow_rate(flow_rate_m3_s: VolumetricFlowRate<f64>) -> Self {
        Self {
            boundary: BranchBoundarySpecification::FlowRate { flow_rate_m3_s },
        }
    }
}

crate::impl_metadata!(BranchBoundaryMetadata, "BranchBoundaryMetadata");

#[cfg(test)]
mod tests {
    use super::{BranchBoundaryMetadata, BranchBoundarySpecification};
    use aequitas::systems::si::quantities::{Pressure, VolumetricFlowRate};

    #[test]
    fn branch_boundary_metadata_preserves_typed_values_through_json() {
        let authored = BranchBoundaryMetadata::pressure(Pressure::from_base(125.0));
        let encoded = serde_json::to_string(&authored).expect("metadata serializes");
        let decoded: BranchBoundaryMetadata =
            serde_json::from_str(&encoded).expect("metadata deserializes");

        assert_eq!(decoded, authored);
        assert_eq!(
            decoded.boundary,
            BranchBoundarySpecification::Pressure {
                pressure_pa: Pressure::from_base(125.0),
            }
        );

        let flow = BranchBoundaryMetadata::flow_rate(VolumetricFlowRate::from_base(-2.5e-7));
        let flow_encoded = serde_json::to_string(&flow).expect("flow serializes");
        let decoded_flow: BranchBoundaryMetadata =
            serde_json::from_str(&flow_encoded).expect("flow deserializes");
        assert_eq!(decoded_flow, flow);
    }
}

/// Optimization metadata recorded for a length-optimized channel.
#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationMetadata {
    /// Authored channel length before optimization.
    pub original_length: f64,
    /// Channel length after optimization.
    pub optimized_length: f64,
    /// Relative length improvement as a percentage.
    pub improvement_percentage: f64,
    /// Number of optimization iterations.
    pub iterations: usize,
    /// Wall-clock optimization time in milliseconds.
    pub optimization_time_ms: u64,
    /// Optimization profile identifier.
    pub optimization_profile: String,
}

crate::impl_metadata!(OptimizationMetadata, "OptimizationMetadata");

/// Performance metadata recorded for blueprint generation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerformanceMetadata {
    /// Generation time in microseconds.
    pub generation_time_us: u64,
    /// Peak memory usage in bytes.
    pub memory_usage_bytes: usize,
    /// Number of path points in the generated geometry.
    pub path_points_count: usize,
}

crate::impl_metadata!(PerformanceMetadata, "PerformanceMetadata");
