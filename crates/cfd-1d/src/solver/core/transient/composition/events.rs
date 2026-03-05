use super::state::MixtureComposition;
use nalgebra::RealField;

/// Piecewise-constant inlet mixture event.
#[derive(Debug, Clone)]
pub struct InletCompositionEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Node index where the inlet composition applies.
    pub node_index: usize,
    /// Mixture after this event.
    pub mixture: MixtureComposition<T>,
}

/// Piecewise-constant edge flow event for transient pump-style control.
#[derive(Debug, Clone)]
pub struct EdgeFlowEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Edge index whose flow is updated.
    pub edge_index: usize,
    /// Flow rate after this event.
    pub flow_rate: T,
}

/// Piecewise-constant node pressure boundary event for transient pressure-pump control.
#[derive(Debug, Clone)]
pub struct PressureBoundaryEvent<T: RealField + Copy> {
    /// Event activation time.
    pub time: T,
    /// Node index whose pressure boundary is updated.
    pub node_index: usize,
    /// Pressure value after this event.
    pub pressure: T,
}
