//! Unified fluid property traits
//!
//! This module provides the single source of truth for fluid property interfaces,
//! reconciling the needs of both constant and variable property models.

use crate::error::Error;
use crate::physics::fluid::thermophysical;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Length, MassDensity, SpecificHeatCapacity,
    ThermalConductivity, Velocity,
};
use eunomia::NumericElement;
use eunomia::RealField;

/// Core fluid properties structure
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FluidState<T: RealField + Copy> {
    /// Density [kg/m³]
    pub density: T,
    /// Dynamic viscosity [Pa·s]
    pub dynamic_viscosity: T,
    /// Specific heat capacity [J/(kg·K)]
    pub specific_heat: T,
    /// Thermal conductivity [W/(m·K)]
    pub thermal_conductivity: T,
    /// Speed of sound \[m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + NumericElement + Copy> FluidState<T> {
    /// Calculate kinematic viscosity [m²/s]
    #[inline]
    pub fn kinematic_viscosity(&self) -> T {
        let dynamic_viscosity = DynamicViscosity::from_base(self.dynamic_viscosity);
        let density = MassDensity::from_base(self.density);
        let kinematic = dynamic_viscosity / density;
        kinematic.into_base()
    }

    /// Calculate Prandtl number
    #[inline]
    pub fn prandtl_number(&self) -> T {
        let dynamic_viscosity = DynamicViscosity::from_base(self.dynamic_viscosity);
        let specific_heat = SpecificHeatCapacity::from_base(self.specific_heat);
        let thermal_conductivity = ThermalConductivity::from_base(self.thermal_conductivity);
        let prandtl: Dimensionless<T> = dynamic_viscosity * specific_heat / thermal_conductivity;
        prandtl.into_base()
    }

    /// Calculate thermal diffusivity [m²/s]
    ///
    /// # Errors
    ///
    /// Returns an error if density, specific heat, or thermal conductivity
    /// violates the Proteus thermophysical-property contract.
    #[inline]
    pub fn thermal_diffusivity(&self) -> Result<T, Error> {
        thermophysical::thermal_diffusivity(
            self.density,
            self.specific_heat,
            self.thermal_conductivity,
        )
    }

    /// Calculate Reynolds number
    #[inline]
    pub fn reynolds_number(&self, velocity: T, length: T) -> T {
        let density = MassDensity::from_base(self.density);
        let velocity = Velocity::from_base(velocity);
        let length = Length::from_base(length);
        let dynamic_viscosity = DynamicViscosity::from_base(self.dynamic_viscosity);
        let reynolds: Dimensionless<T> = density * velocity * length / dynamic_viscosity;
        reynolds.into_base()
    }

    /// Calculate Peclet number
    #[inline]
    pub fn peclet_number(&self, velocity: T, length: T) -> T {
        self.reynolds_number(velocity, length) * self.prandtl_number()
    }

    /// Calculate Mach number
    #[inline]
    pub fn mach_number(&self, velocity: T) -> T {
        let v_abs = if velocity >= <T as NumericElement>::ZERO {
            velocity
        } else {
            -velocity
        };
        v_abs / self.speed_of_sound
    }
}

/// Unified trait for all fluid models
///
/// This trait provides a single interface for both constant and variable property fluids,
/// eliminating the dual-trait hierarchy that violated SSOT.
pub trait Fluid<T: RealField + Copy>: Send + Sync {
    /// Get fluid properties at given conditions
    fn properties_at(&self, temperature: T, pressure: T) -> Result<FluidState<T>, Error>;

    /// Get speed of sound at given conditions \[m/s]
    fn speed_of_sound_at(&self, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties_at(temperature, pressure)
            .map(|s| s.speed_of_sound)
    }

    /// Get descriptive name
    fn name(&self) -> &str;

    /// Check if properties vary with temperature
    fn is_temperature_dependent(&self) -> bool {
        false
    }

    /// Check if properties vary with pressure
    fn is_pressure_dependent(&self) -> bool {
        false
    }

    /// Get reference temperature for property evaluation \[K]
    fn reference_temperature(&self) -> Option<T> {
        None
    }

    /// Get reference pressure for property evaluation \[Pa]
    fn reference_pressure(&self) -> Option<T> {
        None
    }

    /// Get viscosity at specific shear rate [Pa·s]
    ///
    /// For Newtonian fluids, this returns the dynamic viscosity.
    /// For non-Newtonian fluids, this returns the apparent viscosity.
    fn viscosity_at_shear(&self, _shear_rate: T, temperature: T, pressure: T) -> Result<T, Error> {
        self.properties_at(temperature, pressure)
            .map(|s| s.dynamic_viscosity)
    }
}

/// Convenience trait for constant property fluids
///
/// Provides direct field access for fluids that don't vary with conditions
pub trait ConstantFluid<T: RealField + Copy>: Fluid<T> {
    /// Get constant density
    fn density(&self) -> T;

    /// Get constant dynamic viscosity
    fn dynamic_viscosity(&self) -> T;

    /// Get constant specific heat
    fn specific_heat(&self) -> T;

    /// Get constant thermal conductivity  
    fn thermal_conductivity(&self) -> T;

    /// Get constant speed of sound
    fn speed_of_sound(&self) -> T;

    /// Get constant kinematic viscosity
    fn kinematic_viscosity(&self) -> T {
        self.dynamic_viscosity() / self.density()
    }
}

/// Trait for non-Newtonian fluids
pub trait NonNewtonianFluid<T: RealField + Copy>: Fluid<T> {
    /// Get apparent viscosity at given shear rate
    fn apparent_viscosity(&self, shear_rate: T) -> T;

    /// Check if fluid exhibits yield stress
    fn has_yield_stress(&self) -> bool {
        false
    }

    /// Get yield stress if applicable
    fn yield_stress(&self) -> Option<T> {
        None
    }
}

/// Trait for compressible fluids
pub trait CompressibleFluid<T: RealField + Copy>: Fluid<T> {
    /// Calculate density from equation of state
    fn density_from_state(&self, temperature: T, pressure: T) -> T;

    /// Calculate speed of sound
    fn speed_of_sound(&self, temperature: T, pressure: T) -> T;

    /// Get specific heat ratio (gamma)
    fn heat_capacity_ratio(&self) -> T;
}

// Note: Individual types must implement the domains trait explicitly to avoid conflicts

#[cfg(test)]
mod tests {
    use super::FluidState;

    #[test]
    fn typed_fluid_numbers_preserve_dimensionless_oracles() {
        let state = FluidState {
            density: 1_000.0_f64,
            dynamic_viscosity: 0.004,
            specific_heat: 4_000.0,
            thermal_conductivity: 0.6,
            speed_of_sound: 1_500.0,
        };

        assert_eq!(state.kinematic_viscosity().to_bits(), 4.0e-6_f64.to_bits());
        assert_eq!(
            state.prandtl_number().to_bits(),
            (0.004_f64 * 4_000.0 / 0.6).to_bits()
        );
        assert_eq!(
            state.reynolds_number(0.2, 0.01).to_bits(),
            (1_000.0_f64 * 0.2 * 0.01 / 0.004).to_bits()
        );
    }
}
