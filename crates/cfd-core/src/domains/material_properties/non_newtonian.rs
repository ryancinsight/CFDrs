//! Non-Newtonian fluid models

use super::traits::FluidProperties;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

// Named constants for material properties
const SOLID_LIKE_VISCOSITY: f64 = 1e6; // High viscosity for zero shear rate
const YIELD_STRESS_VISCOSITY: f64 = 1e10; // Very high viscosity below yield stress

/// Power-law fluid model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PowerLawFluid<T: RealField + Copy> {
    /// Fluid density
    pub density: T,
    /// Consistency index
    pub consistency_index: T,
    /// Flow behavior index
    pub flow_behavior_index: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
}

impl<T: RealField + Copy> FluidProperties<T> for PowerLawFluid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn dynamic_viscosity(&self) -> T {
        // For power-law fluids, viscosity depends on shear rate
        // This returns the consistency index as a reference value
        self.consistency_index
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }
}

impl<T: RealField + Copy> PowerLawFluid<T> {
    /// Calculate dynamic viscosity at given shear rate
    /// μ = K * γ^(n-1)
    pub fn dynamic_viscosity_at_shear_rate(&self, shear_rate: T) -> T {
        let shear_rate_abs = shear_rate.abs();
        if shear_rate_abs > T::zero() {
            let n_minus_one = self.flow_behavior_index - T::one();
            self.consistency_index * shear_rate_abs.powf(n_minus_one)
        } else {
            // At zero shear rate, return high viscosity for numerical stability
            T::from_f64(SOLID_LIKE_VISCOSITY)
                .unwrap_or(T::from_f64(SOLID_LIKE_VISCOSITY).unwrap_or_else(|| T::one()))
        }
    }
}

/// Bingham plastic fluid model
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BinghamFluid<T: RealField + Copy> {
    /// Fluid density
    pub density: T,
    /// Plastic viscosity (after yielding)
    pub plastic_viscosity: T,
    /// Yield stress
    pub yield_stress: T,
    /// Thermal conductivity
    pub thermal_conductivity: T,
    /// Specific heat capacity
    pub specific_heat: T,
}

impl<T: RealField + Copy> FluidProperties<T> for BinghamFluid<T> {
    fn density(&self) -> T {
        self.density
    }

    fn dynamic_viscosity(&self) -> T {
        // For Bingham fluids, effective viscosity depends on shear stress
        // This returns the plastic viscosity as a reference value
        self.plastic_viscosity
    }

    fn thermal_conductivity(&self) -> T {
        self.thermal_conductivity
    }

    fn specific_heat(&self) -> T {
        self.specific_heat
    }
}

impl<T: RealField + Copy> BinghamFluid<T> {
    /// Calculate dynamic viscosity at given shear stress
    pub fn dynamic_viscosity_at_shear_stress(&self, shear_stress: T) -> T {
        // Bingham model:
        // If τ < τ_y: material behaves as solid (infinite viscosity)
        // If τ ≥ τ_y: μ = μ_p + τ_y/γ
        let shear_stress_abs = shear_stress.abs();
        if shear_stress_abs < self.yield_stress {
            // Below yield stress - solid-like behavior
            T::from_f64(YIELD_STRESS_VISCOSITY)
                .unwrap_or(T::from_f64(YIELD_STRESS_VISCOSITY).unwrap_or_else(|| T::one()))
        } else {
            // Above yield stress - flows with plastic viscosity
            // Effective viscosity includes yield stress contribution
            // μ_eff = μ_p + τ_y/γ where γ = (τ - τ_y)/μ_p
            let shear_rate = (shear_stress_abs - self.yield_stress) / self.plastic_viscosity;
            if shear_rate > T::zero() {
                self.plastic_viscosity + self.yield_stress / shear_rate
            } else {
                self.plastic_viscosity
            }
        }
    }
}
