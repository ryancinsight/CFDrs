//! Detailed Venturi flow analysis and pressure decomposition.
//!
//! Provides the [`VenturiAnalysis`] result struct with comprehensive
//! breakdown of contraction, friction, expansion, and recovery contributions.

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
        let eta_r = T::from_f64(self.expansion_type.recovery_efficiency()).unwrap_or_else(T::one);

        // ΔP_contraction = ½ρV_t²(1 − β⁴) / C_d²  where β⁴ = (A_t/A_i)² = beta_sq·beta_sq
        let dp_contraction =
            half * density * v_throat * v_throat * (one - beta_sq * beta_sq) / (c_d * c_d);
        let dp_friction =
            f * (self.throat_length / self.throat_diameter) * half * density * v_throat * v_throat;
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
