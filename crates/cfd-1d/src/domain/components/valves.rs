//! Valve components for microfluidic networks
//!
//! # Valve Resistance Theorem
//!
//! A microvalve is characterised by a flow coefficient `Cv` [m³/s/√Pa] such that
//! the pressure drop at full opening follows:
//! ```text
//! ΔP = Q² / Cv²   →   k = 1/Cv²   [Pa·s²/m⁶]
//! ```
//! When partially open with fractional opening `f ∈ [0, 1]`, the effective flow
//! coefficient scales as `Cv_eff = Cv · f`, giving:
//! ```text
//! ΔP = Q² / (Cv · f)²   →   k = 1 / (Cv · f)²
//! ```
//! At `f = 0` (closed) the valve is modelled as a very high linear resistance
//! `R_closed = 10¹²` Pa·s/m³ to prevent division by zero and maintain a
//! well-conditioned matrix.
//!
//! The **linear** resistance coefficient `R` is zero for an open valve because the
//! valve introduces only quadratic (dynamic pressure) losses with no viscous dissipation
//! at the device level.
//!
//! # Invariants
//! - `opening ∈ [0, 1]`: clamped on set
//! - `Cv > 0`: physically required (positive conductance)

use super::{real_from_f64, Component};
use cfd_core::error::Result;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Linear resistance applied when a valve is fully closed [Pa·s/m³].
///
/// This large finite value (10¹²) prevents division by zero while keeping
/// the conductance matrix well-conditioned. A closed valve effectively
/// blocks flow: at 1 µL/s the ΔP would be 10⁶ Pa, far exceeding any
/// realistic microfluidic pressure budget.
const CLOSED_VALVE_RESISTANCE: f64 = 1e12;

/// Valve type enumeration
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ValveType {
    /// Normally open valve
    NormallyOpen,
    /// Normally closed valve
    NormallyClosed,
    /// Check valve (one-way)
    Check,
    /// Proportional control valve
    Proportional,
}

/// Microvalve component
///
/// Models a proportional valve with flow coefficient `Cv` and fractional opening.
/// See module docs for the full theorem.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Microvalve<T: RealField + Copy> {
    /// Flow coefficient [m³/s/√Pa] at full opening
    pub cv: T,
    /// Opening fraction (0 = closed, 1 = fully open; clamped to [0, 1])
    pub opening: T,
    /// Additional parameters
    pub parameters: HashMap<String, T>,
}

impl<T: RealField + Copy + FromPrimitive> Microvalve<T> {
    /// Create a new microvalve with specified flow coefficient (fully open)
    pub fn new(cv: T) -> Self {
        Self {
            cv,
            opening: T::one(),
            parameters: HashMap::new(),
        }
    }

    /// Effective flow coefficient at current opening: `Cv_eff = Cv * opening`
    pub fn effective_cv(&self) -> T {
        self.cv * self.opening
    }
}

impl<T: RealField + Copy + FromPrimitive> Component<T> for Microvalve<T> {
    /// Returns the linear hydraulic resistance.
    ///
    /// - **Closed** (`opening ≤ 0`): large linear resistance `1e12` Pa·s/m³ prevents
    ///   division-by-zero and maintains matrix conditioning.
    /// - **Open**: `0.0` — all losses are quadratic (captured by `coefficients`).
    fn resistance(&self, _fluid: &ConstantPropertyFluid<T>) -> T {
        if self.opening <= T::zero() {
            real_from_f64(CLOSED_VALVE_RESISTANCE)
        } else {
            T::zero()
        }
    }

    fn coefficients(&self, fluid: &ConstantPropertyFluid<T>) -> (T, T) {
        if self.opening <= T::zero() {
            (real_from_f64(CLOSED_VALVE_RESISTANCE), T::zero())
        } else {
            // k = 1 / (Cv * opening)^2
            let cv_eff = self.effective_cv();
            let k = T::one() / (cv_eff * cv_eff);
            (self.resistance(fluid), k)
        }
    }

    fn component_type(&self) -> &'static str {
        "Microvalve"
    }

    fn parameters(&self) -> &HashMap<String, T> {
        &self.parameters
    }

    fn set_parameter(&mut self, key: &str, value: T) -> Result<()> {
        match key {
            "cv" => self.cv = value,
            "opening" => {
                self.opening = if value < T::zero() {
                    T::zero()
                } else if value > T::one() {
                    T::one()
                } else {
                    value
                };
            }
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
    fn test_valve_fully_closed_high_resistance() {
        let fluid = water_20c::<f64>().unwrap();
        let mut valve = Microvalve::new(0.01_f64);
        valve.set_parameter("opening", 0.0).unwrap();
        assert!(
            valve.resistance(&fluid) > 1e10,
            "Closed valve must have very high resistance"
        );
    }

    #[test]
    fn test_valve_open_zero_linear_resistance() {
        let fluid = water_20c::<f64>().unwrap();
        let valve = Microvalve::new(0.01_f64);
        assert_relative_eq!(valve.resistance(&fluid), 0.0, epsilon = 1e-30);
    }

    #[test]
    fn test_valve_coefficients_full_open() {
        let fluid = water_20c::<f64>().unwrap();
        let cv = 0.01_f64;
        let valve = Microvalve::new(cv);
        let (r, k) = valve.coefficients(&fluid);
        assert_relative_eq!(r, 0.0, epsilon = 1e-30);
        assert_relative_eq!(k, 1.0 / (cv * cv), epsilon = 1e-15);
    }

    #[test]
    fn test_valve_coefficients_partial_open() {
        let fluid = water_20c::<f64>().unwrap();
        let cv = 0.01_f64;
        let opening = 0.5_f64;
        let mut valve = Microvalve::new(cv);
        valve.set_parameter("opening", opening).unwrap();
        let (_, k) = valve.coefficients(&fluid);
        let expected_k = 1.0 / (cv * opening).powi(2);
        assert_relative_eq!(k, expected_k, epsilon = 1e-15);
    }

    #[test]
    fn test_valve_opening_clamped_above_one() {
        let mut valve = Microvalve::<f64>::new(0.01);
        valve.set_parameter("opening", 1.5).unwrap();
        assert_relative_eq!(valve.opening, 1.0, epsilon = 1e-15);
    }

    #[test]
    fn test_valve_opening_clamped_below_zero() {
        let mut valve = Microvalve::<f64>::new(0.01);
        valve.set_parameter("opening", -0.5).unwrap();
        assert_relative_eq!(valve.opening, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn test_closed_valve_resistance_approx_1e12() {
        let fluid = water_20c::<f64>().unwrap();
        let mut valve = Microvalve::new(0.01_f64);
        valve.set_parameter("opening", 0.0).unwrap();
        let r = valve.resistance(&fluid);
        assert_relative_eq!(r, 1e12, epsilon = 1.0);
    }

    #[test]
    fn test_open_valve_quadratic_coefficient_from_cv() {
        let fluid = water_20c::<f64>().unwrap();
        let cv = 0.02_f64;
        let valve = Microvalve::new(cv);
        let (r, k) = valve.coefficients(&fluid);
        // Open valve: zero linear resistance, quadratic = 1/Cv^2
        assert_relative_eq!(r, 0.0, epsilon = 1e-30);
        assert_relative_eq!(k, 1.0 / (cv * cv), epsilon = 1e-10);
    }

    #[test]
    fn test_partial_opening_intermediate_resistance() {
        let fluid = water_20c::<f64>().unwrap();
        let cv = 0.01_f64;
        let mut valve = Microvalve::new(cv);
        valve.set_parameter("opening", 0.5).unwrap();
        let (_, k_half) = valve.coefficients(&fluid);

        let valve_full = Microvalve::new(cv);
        let (_, k_full) = valve_full.coefficients(&fluid);

        // Half-open should have 4x higher quadratic coefficient than full-open
        // k = 1/(Cv*f)^2, so k(0.5)/k(1.0) = 4
        assert_relative_eq!(k_half / k_full, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_component_type_is_microvalve() {
        let valve = Microvalve::<f64>::new(0.01);
        assert_eq!(valve.component_type(), "Microvalve");
    }
}
