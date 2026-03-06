//! Flow state and channel model definitions
//!
//! ## Theorem: Knudsen Number Regime Classification
//!
//! The Knudsen number `Kn = λ / Dh` classifies rarefaction effects:
//!
//! | Kn range | Regime |
//! |----------|--------|
//! | < 0.001  | Continuum (Stokes/Laminar/Turbulent) |
//! | 0.001–0.1| Slip flow — velocity slip at wall, non-zero tangential velocity |
//! | 0.1–10   | Transition flow |
//! | > 10     | Free molecular flow |
//!
//! **Reference**: Schaaf, S. A. & Chambre, P. L. (1961). *Flow of Rarefied Gases*.
//! Princeton University Press.
//!
//! At millifluidic scales (Dh ≈ 0.1–1 mm), Kn is always well below 0.001 for liquids,
//! so `classify_flow_regime` only reaches `SlipFlow` for gases at very low pressures
//! or when the channel is nanoscale.

use super::geometry::ChannelGeometry;
use nalgebra::RealField;

/// Minimum Knudsen number for slip-flow regime (Schaaf-Chambre 1961)
pub const KN_SLIP_MIN: f64 = 0.001;

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
    /// Knudsen number `Kn = λ / Dh` — `None` if not computed
    ///
    /// Used to classify slip-flow conditions for rarefied gases.
    pub knudsen_number: Option<T>,
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
    /// Stokes flow (Re < 1, creeping flow)
    Stokes,
    /// Laminar flow (1 ≤ Re < 2300)
    Laminar,
    /// Transitional flow (2300 ≤ Re ≤ 4000)
    Transitional,
    /// Turbulent flow (Re > 4000)
    Turbulent,
    /// Slip flow (0.001 ≤ Kn < 0.1) — rarefied gas or nanoscale liquid
    SlipFlow,
}

impl FlowRegime {
    /// Determine flow regime from Reynolds number alone (no Kn data available).
    ///
    /// Use `classify_with_knudsen` when the Knudsen number is known.
    pub fn from_reynolds_number<T: RealField + Copy + num_traits::FromPrimitive>(re: T) -> Self {
        let re_1 = T::one();
        let re_2300 = T::from_f64(2300.0).expect("Mathematical constant conversion compromised");
        let re_4000 = T::from_f64(4000.0).expect("Mathematical constant conversion compromised");

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

    /// Classify regime with Knudsen number taking priority for slip/rarefied conditions.
    ///
    /// `Kn = λ / Dh` where `λ` is the fluid mean free path and `Dh` is the hydraulic diameter.
    /// When `Kn ≥ 0.001`, slip flow overrides the continuum classification.
    pub fn classify_with_knudsen<
        T: RealField + Copy + num_traits::FromPrimitive + num_traits::ToPrimitive,
    >(
        re: T,
        kn: T,
    ) -> Self {
        let kn_val = kn.to_f64().unwrap_or(0.0);
        if kn_val >= KN_SLIP_MIN {
            FlowRegime::SlipFlow
        } else {
            Self::from_reynolds_number(re)
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
