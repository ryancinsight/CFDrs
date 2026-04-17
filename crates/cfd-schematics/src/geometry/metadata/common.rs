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
    pub channel_diameter_mm: f64,
}

crate::impl_metadata!(ChannelGeometryMetadata, "ChannelGeometryMetadata");

/// Branch-level boundary condition metadata authored in schematics.
///
/// This metadata lets a blueprint mark an endpoint branch with an explicit
/// pressure or flow boundary condition instead of relying only on the legacy
/// `NodeKind` inference in reduced-order solvers.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum BranchBoundarySpecification {
    /// Fixed pressure boundary [Pa].
    Pressure { pressure_pa: f64 },
    /// Fixed volumetric flow boundary [m^3/s].
    ///
    /// The sign convention matches the 1D solver: positive values act as a
    /// source term, negative values as a sink term.
    FlowRate { flow_rate_m3_s: f64 },
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
    pub fn pressure(pressure_pa: f64) -> Self {
        Self {
            boundary: BranchBoundarySpecification::Pressure { pressure_pa },
        }
    }

    /// Create a fixed-flow branch boundary.
    #[must_use]
    pub fn flow_rate(flow_rate_m3_s: f64) -> Self {
        Self {
            boundary: BranchBoundarySpecification::FlowRate { flow_rate_m3_s },
        }
    }
}

crate::impl_metadata!(BranchBoundaryMetadata, "BranchBoundaryMetadata");

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
