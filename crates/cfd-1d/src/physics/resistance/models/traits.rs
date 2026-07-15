//! Traits and common types for resistance models.

use crate::scalar::Cfd1dScalar;
use cfd_core::error::Result;
use cfd_core::physics::fluid::FluidTrait;
use eunomia::{FloatElement, NumericElement};

/// Scalar contract for hydraulic resistance models.
///
/// `cfd-core::physics::fluid::FluidTrait` still requires `Cfd1dScalar`; Eunomia
/// owns scalar construction and scalar-to-f64 diagnostics for resistance math.
pub trait ResistanceScalar: Cfd1dScalar + FloatElement + Copy {}

impl<T> ResistanceScalar for T where T: Cfd1dScalar + FloatElement + Copy {}

#[inline]
pub(crate) fn scalar_from_f64<T: ResistanceScalar>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

#[inline]
pub(crate) fn scalar_to_f64<T: NumericElement>(value: T) -> f64 {
    <T as NumericElement>::to_f64(value)
}

/// Trait for hydraulic resistance models
pub trait ResistanceModel<T: ResistanceScalar> {
    /// Calculate hydraulic resistance [Pa·s/m³]
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T>;

    /// Calculate linear and quadratic resistance coefficients
    ///
    /// Returns (R, k) where ΔP = R·Q + k·Q|Q|
    /// - R: Linear resistance coefficient [Pa·s/m³]
    /// - k: Quadratic loss coefficient [Pa·s²/m⁶]
    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
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
    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)
    }

    /// Validate Mach number for incompressibility assumption (Ma < 0.3)
    fn validate_mach_number<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Mach number validation: Ma < 0.3 for incompressibility
        if let Some(velocity) = conditions.velocity {
            let speed_of_sound =
                fluid.speed_of_sound_at(conditions.temperature, conditions.pressure)?;
            if speed_of_sound > T::zero() {
                let v_abs = if velocity >= T::zero() {
                    velocity
                } else {
                    -velocity
                };
                let mach = v_abs / speed_of_sound;
                let mach_limit = scalar_from_f64::<T>(0.3);
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
pub struct FlowConditions<T> {
    /// Reynolds number
    pub reynolds_number: Option<T>,
    /// Flow velocity \[m/s]
    pub velocity: Option<T>,
    /// Flow rate [m³/s]
    pub flow_rate: Option<T>,
    /// Shear rate [1/s] (at wall or characteristic)
    pub shear_rate: Option<T>,
    /// Temperature \[K]
    pub temperature: T,
    /// Pressure \[Pa]
    pub pressure: T,
}

impl<T: ResistanceScalar> FlowConditions<T> {
    /// Create new flow conditions with default temperature and pressure
    pub fn new(velocity: T) -> Self {
        use cfd_core::physics::constants::physics::thermo::{P_ATM, T_STANDARD};

        Self {
            reynolds_number: None,
            velocity: Some(velocity),
            flow_rate: None,
            shear_rate: None,
            temperature: scalar_from_f64::<T>(T_STANDARD),
            pressure: scalar_from_f64::<T>(P_ATM),
        }
    }

    /// Create flow conditions from a prescribed volumetric flow rate.
    ///
    /// The velocity is intentionally left unset so model-specific kinematics can
    /// derive it from geometry rather than being pinned to an explicit zero.
    pub fn from_flow_rate(flow_rate: T) -> Self {
        let mut conditions = Self::new(T::zero());
        conditions.velocity = None;
        conditions.flow_rate = Some(flow_rate);
        conditions
    }
}
