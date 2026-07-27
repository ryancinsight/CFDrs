use super::state::MixtureComposition;
use crate::scalar::Cfd1dScalar;
use aequitas::systems::si::quantities::{Dimensionless, Pressure, Time, VolumetricFlowRate};

/// Piecewise-constant inlet hematocrit event for blood-fraction transport.
#[derive(Debug, Clone)]
pub struct InletHematocritEvent<T: Cfd1dScalar + Copy> {
    /// Event activation time.
    pub time: Time<T>,
    /// Node index where the inlet hematocrit applies.
    pub node_index: usize,
    /// Feed hematocrit after this event.
    pub hematocrit: Dimensionless<T>,
}

/// Piecewise-constant inlet mixture event.
#[derive(Debug, Clone)]
pub struct InletCompositionEvent<T: Cfd1dScalar + Copy> {
    /// Event activation time.
    pub time: Time<T>,
    /// Node index where the inlet composition applies.
    pub node_index: usize,
    /// Mixture after this event.
    pub mixture: MixtureComposition<T>,
}

/// Piecewise-constant edge flow event for transient pump-style control.
#[derive(Debug, Clone)]
pub struct EdgeFlowEvent<T: Cfd1dScalar + Copy> {
    /// Event activation time.
    pub time: Time<T>,
    /// Edge index whose flow is updated.
    pub edge_index: usize,
    /// Flow rate after this event.
    pub flow_rate: VolumetricFlowRate<T>,
}

/// Piecewise-constant node pressure boundary event for transient pressure-pump control.
#[derive(Debug, Clone)]
pub struct PressureBoundaryEvent<T: Cfd1dScalar + Copy> {
    /// Event activation time.
    pub time: Time<T>,
    /// Node index whose pressure boundary is updated.
    pub node_index: usize,
    /// Pressure value after this event.
    pub pressure: Pressure<T>,
}
