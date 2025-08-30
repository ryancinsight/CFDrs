//! Darcy-Weisbach resistance model for turbulent flow.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for power operations
const POWER_TWO: f64 = 2.0;
const AREA_DIVISOR: f64 = 4.0;
const RESISTANCE_DIVISOR: f64 = 2.0;

// Named constants for friction factor calculations
const LAMINAR_FRICTION_COEFFICIENT: f64 = 64.0;
const LAMINAR_TRANSITION_RE: f64 = 2000.0;
const HAALAND_ROUGHNESS_DIVISOR: f64 = 3.6;
const HAALAND_REYNOLDS_FACTOR: f64 = 6.9;
const LOG_BASE_10: f64 = 10.0;
const HAALAND_EXPONENT_FACTOR: f64 = 1.8;
const COLEBROOK_COEFFICIENT: f64 = 2.0;
const MAX_REYNOLDS: f64 = 1e8;

/// Darcy-Weisbach resistance model for turbulent flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DarcyWeisbachModel<T: RealField + Copy> {
    /// Hydraulic diameter [m]
    pub hydraulic_diameter: T,
    /// Channel length [m]
    pub length: T,
    /// Surface roughness [m]
    pub roughness: T,
}

impl<T: RealField + Copy> DarcyWeisbachModel<T> {
    /// Create a new Darcy-Weisbach model
    pub fn new(hydraulic_diameter: T, length: T, roughness: T) -> Self {
        Self {
            hydraulic_diameter,
            length,
            roughness,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for DarcyWeisbachModel<T>
{
    fn calculate_resistance(&self, fluid: &Fluid<T>, conditions: &FlowConditions<T>) -> Result<T> {
        let reynolds = conditions.reynolds_number.ok_or_else(|| {
            Error::InvalidConfiguration(
                "Reynolds number required for Darcy-Weisbach model".to_string(),
            )
        })?;

        // Calculate friction factor using Colebrook-White equation iteratively
        let friction_factor = self.calculate_friction_factor(reynolds);

        let area = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero())
            * num_traits::Float::powf(
                self.hydraulic_diameter,
                T::from_f64(POWER_TWO).unwrap_or_else(|| T::zero()),
            )
            / T::from_f64(AREA_DIVISOR).unwrap_or_else(|| T::zero());

        // Convert Darcy friction factor to hydraulic resistance
        let density = fluid.density;
        let resistance = friction_factor * self.length * density
            / (T::from_f64(RESISTANCE_DIVISOR).unwrap_or_else(|| T::zero())
                * area
                * num_traits::Float::powf(
                    self.hydraulic_diameter,
                    T::from_f64(POWER_TWO).unwrap_or_else(|| T::zero()),
                ));

        Ok(resistance)
    }

    fn model_name(&self) -> &str {
        "Darcy-Weisbach"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::from_f64(cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_UPPER)
                .unwrap_or_else(|| T::zero()),
            T::from_f64(MAX_REYNOLDS).unwrap_or_else(|| T::zero()),
        )
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> DarcyWeisbachModel<T> {
    /// Calculate friction factor using iterative Colebrook-White equation
    fn calculate_friction_factor(&self, reynolds: T) -> T {
        use crate::components::constants::COLEBROOK_TOLERANCE;
        use cfd_core::constants::physics::hydraulics::{
            COLEBROOK_REYNOLDS_NUMERATOR, COLEBROOK_ROUGHNESS_DIVISOR,
        };

        let relative_roughness = self.roughness / self.hydraulic_diameter;

        // Check for laminar flow
        let re_transition = T::from_f64(LAMINAR_TRANSITION_RE).unwrap_or_else(|| T::one());
        if reynolds < re_transition {
            // Laminar flow: f = 64/Re
            return T::from_f64(LAMINAR_FRICTION_COEFFICIENT).unwrap_or_else(|| T::one())
                / reynolds;
        }

        // Initial guess using Haaland explicit formula for convergence
        let mut f = {
            let term = relative_roughness
                / T::from_f64(HAALAND_ROUGHNESS_DIVISOR).unwrap_or_else(|| T::one())
                + T::from_f64(HAALAND_REYNOLDS_FACTOR).unwrap_or_else(|| T::one()) / reynolds;
            let log_term = num_traits::Float::ln(term)
                / T::from_f64(LOG_BASE_10.ln()).unwrap_or_else(|| T::one());
            T::one()
                / num_traits::Float::powi(
                    T::from_f64(HAALAND_EXPONENT_FACTOR).unwrap_or_else(|| T::one()) * log_term,
                    2,
                )
        };

        // Iterative solution of Colebrook-White equation
        // 1/sqrt(f) = -2.0 * log10(ε/(3.7*D) + 2.51/(Re*sqrt(f)))
        let tolerance = T::from_f64(COLEBROOK_TOLERANCE)
            .unwrap_or_else(|| T::from_f64(1e-6).unwrap_or_else(|| T::zero()));
        let max_iter = 50;
        let ln10 = T::from_f64(LOG_BASE_10.ln()).unwrap_or_else(|| T::one());

        for _ in 0..max_iter {
            let sqrt_f = num_traits::Float::sqrt(f);
            let term1 = relative_roughness
                / T::from_f64(COLEBROOK_ROUGHNESS_DIVISOR).unwrap_or_else(|| T::one());
            let term2 = T::from_f64(COLEBROOK_REYNOLDS_NUMERATOR).unwrap_or_else(|| T::one())
                / (reynolds * sqrt_f);

            let log_arg = term1 + term2;
            if log_arg <= T::zero() {
                // Fallback to previous value if log argument invalid
                break;
            }

            let inv_sqrt_f = -T::from_f64(COLEBROOK_COEFFICIENT).unwrap_or_else(|| T::one())
                * (num_traits::Float::ln(log_arg) / ln10);
            let f_next = T::one() / (inv_sqrt_f * inv_sqrt_f);

            // Check convergence
            if num_traits::Float::abs(f_next - f) < tolerance {
                f = f_next;
                break;
            }

            f = f_next;
        }

        f
    }
}
