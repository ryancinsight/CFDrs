//! Detailed Venturi flow analysis and pressure decomposition.
//!
//! Provides the [`VenturiAnalysis`] result struct with comprehensive
//! breakdown of contraction, friction, expansion, and recovery contributions.
//!
//! # Theorem - Venturi Pressure-Decomposition Identity
//!
//! For a fixed venturi geometry and admissible inlet kinematics, the detailed
//! analysis satisfies
//!
//! ```text
//! ΔP_total = ΔP_contraction + ΔP_friction + ΔP_expansion_loss - ΔP_recovery
//! ```
//!
//! with $ΔP_contraction$ from the Bernoulli acceleration corrected by the
//! discharge coefficient, $ΔP_friction$ from Darcy-Weisbach throat losses,
//! $ΔP_expansion_loss$ from the diffuser loss coefficient, and $ΔP_recovery$
//! from the Reynolds-corrected diffuser recovery efficiency.
//!
//! **Proof sketch**: the implementation evaluates each contribution from the
//! same inlet, throat, and outlet velocities implied by continuity. The net
//! pressure drop is then formed by direct superposition of irreversible losses
//! and reversible recovery, so any equivalent specification of the inlet state
//! through either velocity or flow rate must produce the same decomposition.

use super::model::VenturiModel;
use super::traits::FlowConditions;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// Detailed Venturi flow analysis result
#[derive(Debug, Clone)]
pub struct VenturiAnalysis<T: RealField + Copy> {
    /// Throat velocity [m/s]
    pub throat_velocity: T,
    /// Throat Reynolds number
    pub throat_reynolds: T,
    /// Throat wall shear rate [1/s]
    pub throat_shear_rate: T,
    /// Apparent viscosity at throat [Pa·s]
    pub throat_viscosity: T,
    /// Contraction pressure drop [Pa]
    pub dp_contraction: T,
    /// Throat friction pressure drop [Pa]
    pub dp_friction: T,
    /// Expansion loss [Pa]
    pub dp_expansion_loss: T,
    /// Expansion recovery [Pa]
    pub dp_recovery: T,
    /// Net pressure drop [Pa]
    pub dp_total: T,
    /// Discharge coefficient used
    pub discharge_coefficient: T,
    /// Expansion loss coefficient
    pub expansion_loss_coefficient: T,
    /// Darcy friction factor in throat
    pub friction_factor: T,
}

impl<T: RealField + Copy + FromPrimitive> VenturiModel<T> {
    /// Perform detailed Venturi flow analysis
    ///
    /// Returns comprehensive breakdown of pressure contributions
    pub fn analyze<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<VenturiAnalysis<T>> {
        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let density = state.density;

        let a_inlet = self.inlet_area();
        let a_throat = self.throat_area();
        let a_outlet = self.outlet_area();

        let v_inlet = if let Some(v) = conditions.velocity {
            v
        } else if let Some(q) = conditions.flow_rate {
            q / a_inlet
        } else {
            return Err(Error::InvalidConfiguration(
                "Venturi analysis requires velocity or flow_rate".to_string(),
            ));
        };

        let v_throat = v_inlet * a_inlet / a_throat;
        let v_outlet = v_inlet * a_inlet / a_outlet;

        let eight = T::from_f64(8.0).expect("Mathematical constant conversion compromised");
        let half = T::one() / (T::one() + T::one());
        let one = T::one();

        let shear_rate_throat = eight * v_throat / self.throat_diameter;
        let viscosity = fluid.viscosity_at_shear(
            shear_rate_throat,
            conditions.temperature,
            conditions.pressure,
        )?;

        let re_throat = density * v_throat * self.throat_diameter / viscosity;

        let viscosity_inlet = fluid.viscosity_at_shear(
            eight * v_inlet / self.inlet_diameter,
            conditions.temperature,
            conditions.pressure,
        )?;
        let re_inlet = density * v_inlet * self.inlet_diameter / viscosity_inlet;

        let beta_sq = self.beta_squared();
        let c_d = self.effective_discharge_coefficient(re_inlet);
        let f = self.throat_friction_factor(re_throat);
        let k_exp = T::from_f64(self.expansion_type.loss_coefficient()).unwrap_or_else(T::one);
        let eta_r = self.effective_recovery_efficiency(re_throat);

        // ΔP_contraction = ½ρV_t²(1 − β⁴) / C_d²  where β⁴ = (A_t/A_i)² = beta_sq·beta_sq
        let dp_contraction =
            half * density * v_throat * v_throat * (one - beta_sq * beta_sq) / (c_d * c_d);
        let dp_friction = Self::throat_friction_pressure_drop(
            f,
            density,
            self.throat_length,
            self.throat_diameter,
            v_throat,
        );
        let dv = v_throat - v_outlet;
        let dp_expansion_loss = k_exp * half * density * dv * dv;
        let dp_recovery = eta_r * half * density * (v_throat * v_throat - v_outlet * v_outlet);
        let dp_total = dp_contraction + dp_friction + dp_expansion_loss - dp_recovery;

        Ok(VenturiAnalysis {
            throat_velocity: v_throat,
            throat_reynolds: re_throat,
            throat_shear_rate: shear_rate_throat,
            throat_viscosity: viscosity,
            dp_contraction,
            dp_friction,
            dp_expansion_loss,
            dp_recovery,
            dp_total,
            discharge_coefficient: c_d,
            expansion_loss_coefficient: k_exp,
            friction_factor: f,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::database::water_20c;
    use proptest::prelude::*;

    #[test]
    fn analyze_rejects_missing_velocity_and_flow_rate() -> cfd_core::error::Result<()> {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
        let fluid = water_20c::<f64>()?;
        let mut conditions = FlowConditions::new(0.1);
        conditions.velocity = None;
        conditions.flow_rate = None;

        assert!(model.analyze(&fluid, &conditions).is_err());
        Ok(())
    }

    #[test]
    fn analyze_pressure_drop_matches_reported_components() -> cfd_core::error::Result<()> {
        let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
        let fluid = water_20c::<f64>()?;
        let conditions = FlowConditions::new(0.2);

        let analysis = model.analyze(&fluid, &conditions)?;
        let expected_total =
            analysis.dp_contraction + analysis.dp_friction + analysis.dp_expansion_loss
                - analysis.dp_recovery;

        assert_relative_eq!(
            analysis.dp_total,
            expected_total,
            epsilon = expected_total.abs().max(1.0) * 1e-12
        );
        Ok(())
    }

    proptest! {
        #[test]
        fn analyze_is_invariant_under_equivalent_velocity_and_flow_rate(
            inlet_velocity_m_s in 1.0e-4_f64..2.0,
        ) {
            let model = VenturiModel::symmetric(0.01_f64, 0.005, 0.01, 0.05);
            let fluid = water_20c::<f64>().expect("water properties should load");

            let velocity_conditions = FlowConditions::new(inlet_velocity_m_s);
            let flow_rate_conditions =
                FlowConditions::from_flow_rate(inlet_velocity_m_s * model.inlet_area());

            let velocity_analysis = model
                .analyze(&fluid, &velocity_conditions)
                .expect("velocity-driven venturi analysis should succeed");
            let flow_rate_analysis = model
                .analyze(&fluid, &flow_rate_conditions)
                .expect("flow-rate-driven venturi analysis should succeed");

            assert_relative_eq!(velocity_analysis.throat_velocity, flow_rate_analysis.throat_velocity, epsilon = velocity_analysis.throat_velocity.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.throat_reynolds, flow_rate_analysis.throat_reynolds, epsilon = velocity_analysis.throat_reynolds.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.throat_shear_rate, flow_rate_analysis.throat_shear_rate, epsilon = velocity_analysis.throat_shear_rate.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.throat_viscosity, flow_rate_analysis.throat_viscosity, epsilon = velocity_analysis.throat_viscosity.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.dp_contraction, flow_rate_analysis.dp_contraction, epsilon = velocity_analysis.dp_contraction.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.dp_friction, flow_rate_analysis.dp_friction, epsilon = velocity_analysis.dp_friction.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.dp_expansion_loss, flow_rate_analysis.dp_expansion_loss, epsilon = velocity_analysis.dp_expansion_loss.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.dp_recovery, flow_rate_analysis.dp_recovery, epsilon = velocity_analysis.dp_recovery.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.dp_total, flow_rate_analysis.dp_total, epsilon = velocity_analysis.dp_total.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.discharge_coefficient, flow_rate_analysis.discharge_coefficient, epsilon = velocity_analysis.discharge_coefficient.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.expansion_loss_coefficient, flow_rate_analysis.expansion_loss_coefficient, epsilon = velocity_analysis.expansion_loss_coefficient.abs().max(1.0) * 1e-12);
            assert_relative_eq!(velocity_analysis.friction_factor, flow_rate_analysis.friction_factor, epsilon = velocity_analysis.friction_factor.abs().max(1.0) * 1e-12);
        }
    }
}
