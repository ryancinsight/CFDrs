//! Discretization schemes for momentum equations

use nalgebra::RealField;

/// Trait for discretization schemes
pub trait DiscretizationScheme<T: RealField + Copy> {
    /// Compute convective flux
    fn convective_flux(&self, phi_c: T, phi_u: T, phi_d: T, velocity: T) -> T;

    /// Get scheme name
    fn name(&self) -> &str;
}

/// Upwind discretization scheme
pub struct Upwind;

impl<T: RealField + Copy> DiscretizationScheme<T> for Upwind {
    fn convective_flux(&self, phi_c: T, phi_u: T, _phi_d: T, velocity: T) -> T {
        if velocity > T::zero() {
            velocity * phi_u
        } else {
            velocity * phi_c
        }
    }

    fn name(&self) -> &'static str {
        "Upwind"
    }
}

/// Central difference discretization scheme
pub struct CentralDifference;

impl<T: RealField + Copy> DiscretizationScheme<T> for CentralDifference {
    fn convective_flux(&self, phi_c: T, phi_u: T, _phi_d: T, velocity: T) -> T {
        let half = T::one() / (T::one() + T::one());
        velocity * (phi_c + phi_u) * half
    }

    fn name(&self) -> &'static str {
        "CentralDifference"
    }
}
