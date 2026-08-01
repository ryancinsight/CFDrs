//! Turbulence modeling abstractions
use super::fields::FlowField;
use aequitas::systems::si::quantities::{KinematicViscosity, SpecificEnergy};
use eunomia::NumericElement;

/// Turbulence model abstraction following Strategy pattern
pub trait TurbulenceModel<T: NumericElement>: Send + Sync {
    /// Calculate turbulent kinematic viscosity as a field in m²/s.
    fn turbulent_viscosity(&self, flow_field: &FlowField<T>) -> Vec<KinematicViscosity<T>>;

    /// Calculate turbulent kinetic energy as a specific-energy field in J/kg.
    ///
    /// The scalar is extracted only inside formulas or explicit serialization
    /// boundaries; callers cannot accidentally combine this field with a
    /// dimensionally unrelated scalar field.
    fn turbulent_kinetic_energy(&self, flow_field: &FlowField<T>) -> Vec<SpecificEnergy<T>>;

    /// Get model name
    fn name(&self) -> &str;
}
