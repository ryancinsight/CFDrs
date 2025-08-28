//! Entrance effects resistance model.

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::Result;
use cfd_core::fluid::Fluid;
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constant for maximum Reynolds number
const MAX_REYNOLDS_ENTRANCE: f64 = 1e6;

/// Entrance effects resistance model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EntranceEffectsModel<T: RealField + Copy> {
    /// Base resistance value
    pub base_resistance: T,
    /// Entrance length coefficient
    pub entrance_coefficient: T,
}

impl<T: RealField + Copy + FromPrimitive + num_traits::Float> ResistanceModel<T>
    for EntranceEffectsModel<T>
{
    fn calculate_resistance(
        &self,
        _fluid: &Fluid<T>,
        _conditions: &FlowConditions<T>,
    ) -> Result<T> {
        Ok(self.base_resistance * (T::one() + self.entrance_coefficient))
    }

    fn model_name(&self) -> &str {
        "Entrance Effects"
    }

    fn reynolds_range(&self) -> (T, T) {
        (
            T::zero(),
            T::from_f64(MAX_REYNOLDS_ENTRANCE).unwrap_or_else(|| T::zero()),
        )
    }
}

/// Methods for combining multiple resistance models
#[derive(Debug, Clone, PartialEq)]
pub enum CombinationMethod {
    /// Add resistances in series
    Series,
    /// Add resistances in parallel
    Parallel,
    /// Custom weighted combination
    Weighted(Vec<f64>),
}
