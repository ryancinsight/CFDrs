//! Hagen-Poiseuille resistance model for laminar flow in circular pipes.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::fluid::{ConstantPropertyFluid, Fluid};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants
const HAGEN_POISEUILLE_COEFFICIENT: f64 = 128.0;
const POWER_FOUR: f64 = 4.0;

/// Hagen-Poiseuille resistance model for circular channels
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HagenPoiseuilleModel<T: RealField + Copy> {
    /// Channel diameter [m]
    pub diameter: T,
    /// Channel length [m]
    pub length: T,
}

impl<T: RealField + Copy> HagenPoiseuilleModel<T> {
    /// Create a new Hagen-Poiseuille model
    pub fn new(diameter: T, length: T) -> Self {
        Self { diameter, length }
    }
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for HagenPoiseuilleModel<T>
{
    fn calculate_resistance(&self, fluid: &Fluid<T>, _conditions: &FlowConditions<T>) -> Result<T> {
        let viscosity = fluid.dynamic_viscosity();
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(|| T::zero());

        let coefficient = T::from_f64(HAGEN_POISEUILLE_COEFFICIENT).unwrap_or_else(|| T::zero());

        // R = (128 * μ * L) / (π * D^4)
        let d4 = num_traits::Float::powf(
            self.diameter,
            T::from_f64(POWER_FOUR).unwrap_or_else(|| T::zero()),
        );
        let resistance = coefficient * viscosity * self.length / (pi * d4);

        Ok(resistance)
    }

    fn model_name(&self) -> &str {
        "Hagen-Poiseuille"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(cfd_core::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER)
                .unwrap_or_else(|| T::zero()),
        )
    }
}
