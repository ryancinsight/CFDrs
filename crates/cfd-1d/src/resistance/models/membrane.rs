//! Porous membrane resistance model for organ-on-chip and barrier simulations.
//!
//! Uses a parallel-pore approximation with cylindrical pores:
//! - Single pore: R_pore = 8 μ L / (π r^4)
//! - Number of pores: N = φ A / (π r^2)
//! - Total: R_total = R_pore / N = 8 μ L / (φ A r^2)

use super::traits::{FlowConditions, ResistanceModel};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid::FluidTrait;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Porous membrane represented by equivalent parallel cylindrical pores.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MembranePoreModel<T: RealField + Copy> {
    /// Membrane thickness [m]
    pub thickness: T,
    /// Membrane width [m]
    pub width: T,
    /// Membrane height [m]
    pub height: T,
    /// Pore radius [m]
    pub pore_radius: T,
    /// Open-area fraction, in [0, 1]
    pub porosity: T,
}

impl<T: RealField + Copy> MembranePoreModel<T> {
    /// Create a new membrane pore model.
    pub fn new(thickness: T, width: T, height: T, pore_radius: T, porosity: T) -> Self {
        Self {
            thickness,
            width,
            height,
            pore_radius,
            porosity,
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ResistanceModel<T> for MembranePoreModel<T> {
    fn calculate_resistance<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<T> {
        self.validate_invariants(fluid, conditions)?;

        let state = fluid.properties_at(conditions.temperature, conditions.pressure)?;
        let mu = state.dynamic_viscosity;

        let area = self.width * self.height;
        let denom = self.porosity * area * self.pore_radius * self.pore_radius;
        let eight = T::from_f64(8.0).unwrap_or_else(T::zero);
        Ok(eight * mu * self.thickness / denom)
    }

    fn model_name(&self) -> &str {
        "MembranePore"
    }

    fn reynolds_range(&self) -> (T, T) {
        (T::zero(), T::from_f64(100.0).unwrap_or_else(T::one))
    }

    fn validate_invariants<F: FluidTrait<T>>(
        &self,
        fluid: &F,
        conditions: &FlowConditions<T>,
    ) -> Result<()> {
        self.validate_mach_number(fluid, conditions)?;

        let zero = T::zero();
        let one = T::one();
        if self.thickness <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane thickness must be positive".to_string(),
            ));
        }
        if self.width <= zero || self.height <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane width and height must be positive".to_string(),
            ));
        }
        if self.pore_radius <= zero {
            return Err(Error::InvalidConfiguration(
                "Membrane pore radius must be positive".to_string(),
            ));
        }
        if self.porosity <= zero || self.porosity > one {
            return Err(Error::InvalidConfiguration(
                "Membrane porosity must be in (0, 1]".to_string(),
            ));
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn membrane_model_returns_positive_resistance() {
        let fluid = cfd_core::physics::fluid::database::water_20c::<f64>().expect("water");
        let model = MembranePoreModel::new(10e-6, 1e-3, 1e-3, 0.5e-6, 0.2);
        let conditions = FlowConditions::new(0.0);
        let r = model
            .calculate_resistance(&fluid, &conditions)
            .expect("resistance");
        assert!(r > 0.0);
    }
}
