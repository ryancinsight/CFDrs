//! Traits and common types for resistance models.

use aequitas::systems::si::quantities::{
    Dimensionless, Pressure, ReciprocalTime, ThermodynamicTemperature, Velocity, VolumetricFlowRate,
};
use aequitas::systems::si::quantities::{HydraulicResistance, QuadraticHydraulicResistance};
use cfd_core::CfdScalar;
use cfd_core::error::Result;
use cfd_core::physics::fluid::FluidTrait;
use eunomia::FloatElement;

/// Trait for hydraulic resistance models.
///
/// The two coefficients occupy different dimensions. Adding them directly
/// without multiplying the quadratic coefficient by flow is invalid.
///
/// ```compile_fail,E0308
/// use aequitas::systems::si::quantities::{
///     HydraulicResistance, QuadraticHydraulicResistance,
/// };
/// let r = HydraulicResistance::from_base(2.0);
/// let k = QuadraticHydraulicResistance::from_base(3.0);
/// let _ = r + k;
/// ```
pub trait ResistanceModel<T: CfdScalar> {
    /// Calculate hydraulic resistance [Pa·s/m³]
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<HydraulicResistance<T>>;

    /// Calculate linear and quadratic resistance coefficients
    ///
    /// Returns (R, k) where ΔP = R·Q + k·Q|Q|
    /// - R: Linear resistance coefficient [Pa·s/m³]
    /// - k: Quadratic loss coefficient [Pa·s²/m⁶]
    fn calculate_coefficients<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<(HydraulicResistance<T>, QuadraticHydraulicResistance<T>)> {
        // Default implementation: assume model is purely linear
        Ok((
            self.calculate_resistance(fluid, conditions)?,
            QuadraticHydraulicResistance::from_base(T::ZERO),
        ))
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

    /// Validate Mach number for incompressibility assumption (Ma ≤ 0.3)
    fn validate_mach_number<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        // Mach number validation: Ma ≤ 0.3 for incompressibility
        if let Some(velocity) = conditions.velocity {
            if !velocity.as_base().is_finite() {
                return Err(cfd_core::error::Error::PhysicsViolation(
                    "Mach number requires a finite flow velocity".into(),
                ));
            }
            let speed_of_sound = fluid.speed_of_sound_at(
                conditions.temperature.into_base(),
                conditions.pressure.into_base(),
            )?;
            if speed_of_sound.into_base() > T::ZERO {
                let mach = velocity.abs() / speed_of_sound;
                let mach_limit = Dimensionless::from_base(<T as FloatElement>::from_f64(0.3));
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
            re.into_base() >= re_min && re.into_base() <= re_max
        } else {
            // Cannot determine applicability without Reynolds number
            false
        }
    }
}

/// Flow conditions for resistance calculations.
///
/// ```compile_fail,E0308
/// use aequitas::systems::si::quantities::VolumetricFlowRate;
/// use cfd_1d::FlowConditions;
/// let _ = FlowConditions::new(VolumetricFlowRate::from_base(0.01));
/// ```
#[derive(Debug, Clone)]
pub struct FlowConditions<T> {
    /// Reynolds number
    pub reynolds_number: Option<Dimensionless<T>>,
    /// Flow velocity \[m/s]
    pub velocity: Option<Velocity<T>>,
    /// Flow rate [m³/s]
    pub flow_rate: Option<VolumetricFlowRate<T>>,
    /// Shear rate [1/s] (at wall or characteristic)
    pub shear_rate: Option<ReciprocalTime<T>>,
    /// Temperature \[K]
    pub temperature: ThermodynamicTemperature<T>,
    /// Absolute pressure for material property queries.
    pub pressure: Pressure<T>,
}

impl<T: CfdScalar> FlowConditions<T> {
    /// Create new flow conditions with default temperature and pressure
    pub fn new(velocity: Velocity<T>) -> Self {
        use cfd_core::physics::constants::physics::thermo::{P_ATM, T_STANDARD};

        Self {
            reynolds_number: None,
            velocity: Some(velocity),
            flow_rate: None,
            shear_rate: None,
            temperature: ThermodynamicTemperature::from_base(<T as FloatElement>::from_f64(
                T_STANDARD,
            )),
            pressure: Pressure::from_base(<T as FloatElement>::from_f64(P_ATM)),
        }
    }

    /// Create flow conditions from a prescribed volumetric flow rate.
    ///
    /// The velocity is intentionally left unset so model-specific kinematics can
    /// derive it from geometry rather than being pinned to an explicit zero.
    pub fn from_flow_rate(flow_rate: VolumetricFlowRate<T>) -> Self {
        let mut conditions = Self::new(Velocity::from_base(T::ZERO));
        conditions.velocity = None;
        conditions.flow_rate = Some(flow_rate);
        conditions
    }
}

#[cfg(test)]
mod tests;
