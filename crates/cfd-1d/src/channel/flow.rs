//! Flow state and channel model definitions

use super::geometry::ChannelGeometry;
use nalgebra::RealField;
/// Extended channel flow model
pub struct Channel<T: RealField + Copy> {
    /// Channel geometry
    pub geometry: ChannelGeometry<T>,
    /// Flow state
    pub flow_state: FlowState<T>,
    /// Numerical parameters
    pub numerical_params: NumericalParameters<T>,
}
/// Flow state information
#[derive(Debug, Clone)]
pub struct FlowState<T: RealField + Copy> {
    /// Reynolds number
    pub reynolds_number: Option<T>,
    /// Flow regime
    pub flow_regime: FlowRegime,
    /// Entrance length effects
    pub entrance_effects: bool,
    /// Secondary flow effects
    pub secondary_flows: bool,
/// Flow regime classification
}

#[derive(Debug, Clone, PartialEq)]
pub enum FlowRegime {
    /// Stokes flow (Re << 1)
    Stokes,
    /// Laminar flow
    Laminar,
    /// Transitional flow
    Transitional,
    /// Turbulent flow
    Turbulent,
    /// Slip flow (rarefied gas)
    SlipFlow,
/// Numerical parameters for advanced modeling
}

pub struct NumericalParameters<T: RealField + Copy> {
    /// Number of discretization points
    pub discretization_points: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Include entrance effects
    /// Include surface tension effects
    pub surface_tension_effects: bool,


}
