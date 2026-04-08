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