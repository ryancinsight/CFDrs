//! Traits and common types for resistance models.

use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Trait for hydraulic resistance models
pub trait ResistanceModel<T: RealField + Copy> {
    /// Calculate hydraulic resistance [Pa·s/m³]
    fn calculate_resistance(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<T>;

    /// Calculate linear and quadratic resistance coefficients
    ///
    /// Returns (R, k) where ΔP = R·Q + k·Q|Q|
    /// - R: Linear resistance coefficient [Pa·s/m³]
    /// - k: Quadratic loss coefficient [Pa·s²/m⁶]
    fn calculate_coefficients(
        &self,
        fluid: &Fluid<T>,
        conditions: &FlowConditions<T>,
    ) -> Result<(T, T)> {
        // Default implementation: assume model is purely linear
        Ok((self.calculate_resistance(fluid, conditions)?, T::zero()))
    }

    /// Get model name
    fn model_name(&self) -> &str;

    /// Get applicable Reynolds number range
    fn reynolds_range(&self) -> (T, T);

    /// Validate physical invariants (e.g., Mach number, entrance length)
    ///
    /// Returns Ok(()) if all invariants are satisfied, or an error describing the violation.
    fn validate_invariants(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<()> {
        self.validate_mach_number(fluid, conditions)
    }

    /// Validate Mach number for incompressibility assumption (Ma < 0.3)
    fn validate_mach_number(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<()> {
        // Mach number validation: Ma < 0.3 for incompressibility
        if let Some(velocity) = conditions.velocity {
            let speed_of_sound = fluid.speed_of_sound;
            if speed_of_sound > T::zero() {
                let v_abs = if velocity >= T::zero() {
                    velocity
                } else {
                    -velocity
                };
                let mach = v_abs / speed_of_sound;
                let mach_limit = T::from_f64(0.3).unwrap_or_else(T::zero);
                if mach > mach_limit {
                    return Err(cfd_core::error::Error::PhysicsViolation(format!(
                        "Mach number violation: Ma > 0.3. Incompressibility assumption invalid for model '{}'",
                        self.model_name()
                    )));
                }
            }
        }
        Ok(())
    }

    /// Check if model is applicable for given conditions
    fn is_applicable(&self, conditions: &FlowConditions<T>) -> bool {
        let (re_min, re_max) = self.reynolds_range();
        if let Some(re) = &conditions.reynolds_number {
            *re >= re_min && *re <= re_max
        } else {
            // Cannot determine applicability without Reynolds number
            false
        }
    }
}

/// Flow conditions for resistance calculations
#[derive(Debug, Clone)]
pub struct FlowConditions<T: RealField + Copy> {
    /// Reynolds number
    pub reynolds_number: Option<T>,
    /// Flow velocity [m/s]
    pub velocity: Option<T>,
    /// Flow rate [m³/s]
    pub flow_rate: Option<T>,
    /// Temperature [K]
    pub temperature: T,
    /// Pressure [Pa]
    pub pressure: T,
}

impl<T: RealField + Copy + FromPrimitive> FlowConditions<T> {
    /// Create new flow conditions with default temperature and pressure
    pub fn new(velocity: T) -> Self {
        use cfd_core::constants::physics::thermo::{P_ATM, T_STANDARD};

        Self {
            reynolds_number: None,
            velocity: Some(velocity),
            flow_rate: None,
            temperature: T::from_f64(T_STANDARD).unwrap_or_else(T::one),
            pressure: T::from_f64(P_ATM).unwrap_or_else(T::one),
        }
    }
}
