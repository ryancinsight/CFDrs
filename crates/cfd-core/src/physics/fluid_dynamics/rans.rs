//! Reynolds-Averaged Navier-Stokes (RANS) turbulence models
use super::fields::FlowField;
use super::turbulence::TurbulenceModel;
use nalgebra::RealField;

/// Base trait for RANS models
pub trait RANSModel<T: RealField + Copy>: TurbulenceModel<T> {
    /// Calculate turbulent dissipation rate
    fn dissipation_rate(&self, flow_field: &FlowField<T>) -> Vec<T>;

    /// Get model constants
    fn constants(&self) -> &dyn std::any::Any;
}
