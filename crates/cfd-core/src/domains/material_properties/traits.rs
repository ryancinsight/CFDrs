//! Core traits for material properties

use nalgebra::RealField;

/// Fluid properties abstraction
pub trait FluidProperties<T: RealField + Copy>: Send + Sync {
    /// Get density
    fn density(&self) -> T;

    /// Get dynamic viscosity
    fn dynamic_viscosity(&self) -> T;

    /// Get kinematic viscosity
    fn kinematic_viscosity(&self) -> T {
        self.dynamic_viscosity() / self.density()
    }

    /// Get thermal conductivity
    fn thermal_conductivity(&self) -> T;

    /// Get specific heat capacity
    fn specific_heat(&self) -> T;

    /// Get thermal diffusivity
    fn thermal_diffusivity(&self) -> T {
        self.thermal_conductivity() / (self.density() * self.specific_heat())
    }

    /// Get Prandtl number
    fn prandtl_number(&self) -> T {
        self.kinematic_viscosity() / self.thermal_diffusivity()
    }
}

/// Solid properties abstraction
pub trait SolidProperties<T: RealField + Copy>: Send + Sync {
    /// Get density
    fn density(&self) -> T;

    /// Get Young's modulus
    fn youngs_modulus(&self) -> T;

    /// Get Poisson's ratio
    fn poissons_ratio(&self) -> T;

    /// Get shear modulus
    fn shear_modulus(&self) -> T {
        let two = T::one() + T::one();
        self.youngs_modulus() / (two * (T::one() + self.poissons_ratio()))
    }

    /// Get thermal conductivity
    fn thermal_conductivity(&self) -> T;

    /// Get specific heat capacity
    fn specific_heat(&self) -> T;

    /// Get thermal expansion coefficient
    fn thermal_expansion(&self) -> T;
}

/// Interface properties abstraction
pub trait InterfaceProperties<T: RealField + Copy>: Send + Sync {
    /// Get surface tension
    fn surface_tension(&self) -> T;

    /// Get contact angle
    fn contact_angle(&self) -> T;

    /// Get adhesion energy
    fn adhesion_energy(&self) -> T {
        self.surface_tension() * self.contact_angle().cos()
    }
}
