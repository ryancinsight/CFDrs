//! Turbulence modeling abstractions
use super::fields::FlowField;
use nalgebra::RealField;

/// Turbulence model abstraction following Strategy pattern
pub trait TurbulenceModel<T: RealField + Copy>: Send + Sync {
    /// Calculate turbulent viscosity
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<T>;

    /// Calculate turbulent kinetic energy
    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<T>;

    /// Get model name
    fn name(&self) -> &str;
}
