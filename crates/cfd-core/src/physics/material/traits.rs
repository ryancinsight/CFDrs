//! Material property traits for solids and interfaces
//!
//! This module provides abstractions for non-fluid material properties
//! used in Fluid-Structure Interaction (FSI) and multi-phase simulations.

use eunomia::{FloatElement, NumericElement};

/// Solid properties abstraction for structural and thermal analysis
pub trait SolidProperties<T: FloatElement + Copy>: Send + Sync {
    /// Get density [kg/m³]
    fn density(&self) -> T;

    /// Get Young's modulus \[Pa]
    fn youngs_modulus(&self) -> T;

    /// Get Poisson's ratio [-]
    fn poissons_ratio(&self) -> T;

    /// Get shear modulus \[Pa]
    fn shear_modulus(&self) -> T {
        let two = <T as NumericElement>::ONE + <T as NumericElement>::ONE;
        self.youngs_modulus() / (two * (<T as NumericElement>::ONE + self.poissons_ratio()))
    }

    /// Get thermal conductivity [W/(m·K)]
    fn thermal_conductivity(&self) -> T;

    /// Get specific heat capacity [J/(kg·K)]
    fn specific_heat(&self) -> T;

    /// Get thermal expansion coefficient [1/K]
    fn thermal_expansion(&self) -> T;
}

/// Interface properties for multi-phase and surface physics
pub trait InterfaceProperties<T: FloatElement + Copy>: Send + Sync {
    /// Get surface tension [N/m]
    fn surface_tension(&self) -> T;

    /// Get contact angle \[rad]
    fn contact_angle(&self) -> T;

    /// Get adhesion energy [J/m²]
    fn adhesion_energy(&self) -> T {
        self.surface_tension() * <T as FloatElement>::cos(self.contact_angle())
    }
}
