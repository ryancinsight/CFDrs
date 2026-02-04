//! Unified fluid property traits
//!
//! This module provides the single source of truth for fluid property interfaces,
//! reconciling the needs of both constant and variable property models.

use crate::error::Error;
use nalgebra::RealField;

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
    /// Speed of sound [m/s]
    pub speed_of_sound: T,
}

impl<T: RealField + Copy> FluidState<T> {
    /// Calculate kinematic viscosity [m²/s]
    #[inline]
    pub fn kinematic_viscosity(&self) -> T {
        self.dynamic_viscosity / self.density
    }

    /// Calculate Prandtl number
    #[inline]
    pub fn prandtl_number(&self) -> T {
        self.dynamic_viscosity * self.specific_heat / self.thermal_conductivity
    }

    /// Calculate thermal diffusivity [m²/s]
    #[inline]
    pub fn thermal_diffusivity(&self) -> T {
        self.thermal_conductivity / (self.density * self.specific_heat)
    }

    /// Calculate Reynolds number
    #[inline]
    pub fn reynolds_number(&self, velocity: T, length: T) -> T {
        self.density * velocity * length / self.dynamic_viscosity
    }

    /// Calculate Peclet number
    #[inline]
    pub fn peclet_number(&self, velocity: T, length: T) -> T {
        self.reynolds_number(velocity, length) * self.prandtl_number()
    }

    /// Calculate Mach number
    #[inline]
    pub fn mach_number(&self, velocity: T) -> T {
        let v_abs = if velocity >= T::zero() {
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

    /// Get speed of sound at given conditions [m/s]
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

    /// Get reference temperature for property evaluation [K]
    fn reference_temperature(&self) -> Option<T> {
        None
    }

    /// Get reference pressure for property evaluation [Pa]
    fn reference_pressure(&self) -> Option<T> {
        None
    }

    /// Get viscosity at specific shear rate [Pa·s]
    ///
    /// For Newtonian fluids, this returns the dynamic viscosity.
    /// For non-Newtonian fluids, this returns the apparent viscosity.
    fn viscosity_at_shear(
        &self,
        _shear_rate: T,
        temperature: T,
        pressure: T,
    ) -> Result<T, Error> {
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
