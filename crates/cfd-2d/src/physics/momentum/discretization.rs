//! Discretization schemes for momentum equations

// MusclReconstruction and MusclOrder integrated into TVD module
use super::muscl::MusclReconstruction;
use nalgebra::RealField;

/// Trait for discretization schemes
pub trait DiscretizationScheme<T: RealField + Copy> {
    /// Compute convective flux
    fn convective_flux(&self, phi_c: T, phi_u: T, phi_d: T, velocity: T) -> T;

    /// Get scheme name
    fn name(&self) -> &str;
}

/// MUSCL-based discretization scheme
pub struct MusclDiscretization<T, M>
where
    T: RealField + Copy,
    M: MusclReconstruction<T>,
{
    muscl: M,
    _phantom: std::marker::PhantomData<T>,
}

impl<T, M> MusclDiscretization<T, M>
where
    T: RealField + Copy,
    M: MusclReconstruction<T>,
{
    /// Create new MUSCL discretization scheme
    pub fn new(muscl: M) -> Self {
        Self {
            muscl,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T, M> DiscretizationScheme<T> for MusclDiscretization<T, M>
where
    T: RealField + Copy,
    M: MusclReconstruction<T>,
{
    fn convective_flux(&self, phi_c: T, phi_u: T, phi_d: T, velocity: T) -> T {
        // For MUSCL, we need interface values reconstructed from cell-centered values
        // TODO: Implement full MUSCL stencil handling for interface reconstruction.

        if velocity > T::zero() {
            // Flow from left to right: use left interface value of upstream cell
            // For positive velocity, upstream is left cell (phi_u)
            // We need to reconstruct the right interface of the upstream cell
            // For simplicity, use MUSCL reconstruction
            let phi_interface = self.muscl.reconstruct_right(phi_c, phi_u, phi_d, None);
            velocity * phi_interface
        } else {
            // Flow from right to left: use right interface value of upstream cell
            // For negative velocity, upstream is right cell (phi_d)
            // We need to reconstruct the left interface of the upstream cell
            let phi_interface = self.muscl.reconstruct_left(phi_u, phi_d, phi_c, None);
            velocity * phi_interface
        }
    }

    fn name(&self) -> &str {
        self.muscl.name()
    }
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
