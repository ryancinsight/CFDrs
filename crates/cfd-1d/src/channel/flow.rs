//! Flow state and channel model definitions

use super::geometry::ChannelGeometry;
use nalgebra::RealField;

/// Extended channel flow model
#[derive(Debug, Clone)]
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
}

/// Flow regime classification
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
}

impl FlowRegime {
    /// Determine flow regime from Reynolds number
    pub fn from_reynolds_number<T: RealField + Copy + num_traits::FromPrimitive>(re: T) -> Self {
        let re_1 = T::from_f64(1.0).unwrap_or_else(T::one);
        let re_2300 =
            T::from_f64(2300.0).unwrap_or_else(|| T::from_usize(2300).unwrap_or_else(T::one));
        let re_4000 =
            T::from_f64(4000.0).unwrap_or_else(|| T::from_usize(4000).unwrap_or_else(T::one));

        if re < re_1 {
            FlowRegime::Stokes
        } else if re < re_2300 {
            FlowRegime::Laminar
        } else if re < re_4000 {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }
}

/// Numerical parameters for advanced modeling
#[derive(Debug, Clone)]
pub struct NumericalParameters<T: RealField + Copy> {
    /// Number of discretization points
    pub discretization_points: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Include entrance effects
    pub entrance_effects: bool,
    /// Include surface tension effects
    pub surface_tension_effects: bool,
}
