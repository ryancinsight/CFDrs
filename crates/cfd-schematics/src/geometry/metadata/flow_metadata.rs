use aequitas::systems::si::quantities::{Length, Pressure, VolumetricFlowRate};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq)]
pub struct FlowMetadata {
    pub flow_rate: f64,
    pub pressure_drop: f64,
    pub reynolds_number: f64,
    pub velocity: f64,
}

crate::impl_metadata!(FlowMetadata, "FlowMetadata");

#[derive(Debug, Clone, PartialEq)]
pub struct ThermalMetadata {
    pub temperature: f64,
    pub heat_transfer_coefficient: f64,
    pub thermal_conductivity: f64,
}

crate::impl_metadata!(ThermalMetadata, "ThermalMetadata");

#[derive(Debug, Clone, PartialEq)]
pub struct ManufacturingMetadata {
    pub width_tolerance: f64,
    pub height_tolerance: f64,
    pub surface_roughness: f64,
    pub manufacturing_method: String,
}

crate::impl_metadata!(ManufacturingMetadata, "ManufacturingMetadata");

#[derive(Debug, Clone, PartialEq)]
pub struct ChannelGeometryMetadata {
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
    Pressure { pressure_pa: Pressure<f64> },
    /// Fixed volumetric flow boundary.
    ///
    /// The sign convention matches the 1D solver: positive values act as a
    /// source term, negative values as a sink term.
    FlowRate {
        flow_rate_m3_s: VolumetricFlowRate<f64>,
    },
}

/// Node-level branch boundary metadata.
///
/// Attach this to a `NodeSpec` in the schematics blueprint when a branch end
/// must override the default inlet/outlet pressure inference.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct BranchBoundaryMetadata {
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

#[derive(Debug, Clone, PartialEq)]
pub struct OptimizationMetadata {
    pub original_length: f64,
    pub optimized_length: f64,
    pub improvement_percentage: f64,
    pub iterations: usize,
    pub optimization_time_ms: u64,
    pub optimization_profile: String,
}

crate::impl_metadata!(OptimizationMetadata, "OptimizationMetadata");

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PerformanceMetadata {
    pub generation_time_us: u64,
    pub memory_usage_bytes: usize,
    pub path_points_count: usize,
}

crate::impl_metadata!(PerformanceMetadata, "PerformanceMetadata");
