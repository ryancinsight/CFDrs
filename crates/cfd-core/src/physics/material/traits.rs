//! Material property traits for solids and interfaces
//!
//! This module provides abstractions for non-fluid material properties
//! used in Fluid-Structure Interaction (FSI) and multi-phase simulations.

use aequitas::systems::si::quantities::{
    Angle, Dimensionless, EnergyPerArea, MassDensity, Pressure, ReciprocalTemperature,
    SpecificHeatCapacity, SurfaceTension, ThermalConductivity,
};
use eunomia::FloatElement;

/// Solid properties abstraction for structural and thermal analysis
pub trait SolidProperties<T: FloatElement + Copy>: Send + Sync {
    /// Get density [kg/m³]
    fn density(&self) -> MassDensity<T>;

    /// Get Young's modulus \[Pa]
    fn youngs_modulus(&self) -> Pressure<T>;

    /// Get Poisson's ratio [-]
    fn poissons_ratio(&self) -> Dimensionless<T>;

    /// Get shear modulus \[Pa]
    ///
    /// Required rather than defaulted: `mu = E / (2 (1 + nu))` is owned by
    /// `proteus::elastic::IsotropicModuli`, and a default here would be a
    /// second copy of that identity. Implementors delegate to the provider.
    fn shear_modulus(&self) -> Pressure<T>;

    /// Get thermal conductivity [W/(m·K)]
    fn thermal_conductivity(&self) -> ThermalConductivity<T>;

    /// Get specific heat capacity [J/(kg·K)]
    fn specific_heat(&self) -> SpecificHeatCapacity<T>;

    /// Get thermal expansion coefficient [1/K]
    fn thermal_expansion(&self) -> ReciprocalTemperature<T>;
}

/// Interface properties for multi-phase and surface physics
pub trait InterfaceProperties<T: FloatElement + Copy>: Send + Sync {
    /// Get surface tension [N/m]
    fn surface_tension(&self) -> SurfaceTension<T>;

    /// Get contact angle \[rad]
    fn contact_angle(&self) -> Angle<T>;

    /// Get adhesion energy [J/m²]
    fn adhesion_energy(&self) -> EnergyPerArea<T> {
        let surface_tension = self.surface_tension().into_base();
        let contact_angle = self.contact_angle().into_base();
        EnergyPerArea::from_base(surface_tension * <T as FloatElement>::cos(contact_angle))
    }
}
