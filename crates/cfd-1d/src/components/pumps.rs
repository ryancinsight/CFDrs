//! Pump components for microfluidic networks
//!
//! # Pump Modelling Theorem
//!
//! A micropump is a **pressure source**, not a passive hydraulic resistor. In the
//! 1D network formulation the system matrix `A x = b` is a Laplacian of conductances.
//! All off-diagonal entries must be non-positive (G_ij ≥ 0) and the diagonal must be
//! strictly positive. A negative resistance would violate positive-definiteness and
//! corrupt the solution.
//!
//! The correct way to model a pump in the conductance-matrix framework is as a
//! **Neumann boundary condition** (prescribed flow source) at the driven node:
//! ```text
//! b_source_node += Q_pump   [m³/s]
//! ```
//! or as a **Dirichlet pressure constraint** at the outlet node.
//!
//! The `Component::resistance()` method on `Micropump` therefore returns **zero** —
//! the pump contributes no passive resistance to the matrix. Its drive is applied
//! separately as a Neumann source term in the right-hand-side vector `b`.
//!
//! # Linear Pump Curve
//!
//! Real pumps follow a pump curve `Q(ΔP)`. For a linear pump curve:
//! ```text
//! Q = Q_max * (1 − ΔP / ΔP_max)
//! ```
//! At the stall point `ΔP = ΔP_max`, Q = 0. At free delivery `ΔP = 0`, Q = Q_max.
//! The hydraulic power delivered is `P_hyd = ΔP · Q`.
//! Pump efficiency: `η = P_hyd / P_input`, typical range 0.1–0.7.

use super::{constants, Component};
use cfd_core::error::Result;
use cfd_core::physics::fluid::Fluid;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Pump type enumeration
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum PumpType {
    /// Syringe pump
    Syringe,
    /// Peristaltic pump
    Peristaltic,
    /// Diaphragm pump
    Diaphragm,
    /// Electroosmotic pump
    Electroosmotic,
}

/// Micropump component
///
/// Models a pump as a **pressure source**. In the conductance-matrix formulation the
/// pump contributes zero passive resistance; its drive is applied as a Neumann source
/// at the driven node. See module docs for the full theorem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Micropump<T: RealField + Copy> {
    /// Maximum flow rate [m³/s] (free-delivery, ΔP = 0)
    pub max_flow_rate: T,
    /// Maximum pressure [Pa] (stall point, Q = 0)
    pub max_pressure: T,
    /// Pump efficiency η [-] (P_hyd / P_input)
    pub efficiency: T,
    /// Operating point on the pump curve (0 = stall, 1 = free delivery)
    pub operating_point: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive + Float> Micropump<T> {
    /// Create a new micropump
    pub fn new(max_flow_rate: T, max_pressure: T) -> Self {
        Self {
            max_flow_rate,
            max_pressure,
            efficiency: T::from_f64(constants::DEFAULT_PUMP_EFFICIENCY).unwrap_or_else(T::one),
            operating_point: T::one(),
            parameters: HashMap::new(),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + Float> Component<T> for Micropump<T> {
    /// Returns **zero** — pumps contribute no passive resistance to the conductance matrix.
    ///
    /// **Invariant**: A pump is a Neumann source (Q_pump injected at the driven node), not
    /// a negative-resistance element. Negative resistance corrupts the matrix Laplacian.
    fn resistance(&self, _fluid: &Fluid<T>) -> T {
        T::zero()
    }

    fn component_type(&self) -> &'static str {
        "Micropump"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "max_flow_rate" => self.max_flow_rate = value,
            "max_pressure" => self.max_pressure = value,
            "efficiency" => self.efficiency = value,
            "operating_point" => self.operating_point = value,
            _ => {
                self.parameters.insert(key.to_string(), value);
            }
        }
        Ok(())
    }

    fn is_active(&self) -> bool {
        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::database::water_20c;

    #[test]
    fn test_pump_resistance_is_zero() {
        let fluid = water_20c::<f64>().unwrap();
        let pump = Micropump::new(1e-6_f64, 1000.0);
        assert_relative_eq!(pump.resistance(&fluid), 0.0, epsilon = 1e-30);
    }

    #[test]
    fn test_pump_is_active() {
        let pump = Micropump::new(1e-6_f64, 1000.0);
        assert!(pump.is_active());
    }

    #[test]
    fn test_pump_parameters_stored() {
        let pump = Micropump::new(1e-6_f64, 1000.0);
        assert_relative_eq!(pump.max_flow_rate, 1e-6, epsilon = 1e-30);
        assert_relative_eq!(pump.max_pressure, 1000.0, epsilon = 1e-30);
    }
}
