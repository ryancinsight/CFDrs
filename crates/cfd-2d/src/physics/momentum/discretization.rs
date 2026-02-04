//! Discretization schemes for momentum equations

// MusclReconstruction and MusclOrder integrated into TVD module
use super::muscl::MusclReconstruction;
use nalgebra::RealField;

/// Trait for discretization schemes
pub trait DiscretizationScheme<T: RealField + Copy> {
    /// Compute convective flux using a 4-point stencil (i-1, i, i+1, i+2)
    /// Flux is computed at the interface between cell i and i+1 (face i+1/2)
    fn convective_flux(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: T, velocity: T) -> T;

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
    fn convective_flux(&self, phi_im1: T, phi_i: T, phi_ip1: T, phi_ip2: T, velocity: T) -> T {
        // For MUSCL, we need interface values reconstructed from cell-centered values
        // Stencil: i-1, i, i+1, i+2
        // Interface: i+1/2 (between i and i+1)

        if velocity > T::zero() {
            // Flow from left to right (positive velocity)
            // We need the Left state of the interface i+1/2 (φ_{i+1/2}^L)
            // This is reconstructed from cell i using upstream (i-1) and downstream (i+1, i+2) neighbors
            // Uses full 4-point stencil: i-1, i, i+1, i+2
            let phi_interface = self
                .muscl
                .reconstruct_left(phi_im1, phi_i, phi_ip1, Some(phi_ip2));
            velocity * phi_interface
        } else {
            // Flow from right to left (negative velocity)
            // We need the Right state of the interface i+1/2 (φ_{i+1/2}^R)
            // This is reconstructed from cell i+1 using upstream (i+2) and downstream (i, i-1) neighbors
            // (Relative to flow direction R->L)
            // For reconstruct_right on cell i+1 (to get left face value), we feed stencil around i+1.
            // Stencil relative to J=i+1: J-1=i, J=i+1, J+1=i+2, J+2=i+3.
            // We have i, i+1, i+2. We are missing i+3.
            // So we pass None for the 4th point.
            let phi_interface = self.muscl.reconstruct_right(phi_i, phi_ip1, phi_ip2, None);
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
    fn convective_flux(&self, _phi_im1: T, phi_i: T, phi_ip1: T, _phi_ip2: T, velocity: T) -> T {
        if velocity > T::zero() {
            // Flow L->R, upwind is i
            velocity * phi_i
        } else {
            // Flow R->L, upwind is i+1
            velocity * phi_ip1
        }
    }

    fn name(&self) -> &'static str {
        "Upwind"
    }
}

/// Central difference discretization scheme
pub struct CentralDifference;

impl<T: RealField + Copy> DiscretizationScheme<T> for CentralDifference {
    fn convective_flux(&self, _phi_im1: T, phi_i: T, phi_ip1: T, _phi_ip2: T, velocity: T) -> T {
        let half = T::one() / (T::one() + T::one());
        velocity * (phi_i + phi_ip1) * half
    }

    fn name(&self) -> &'static str {
        "CentralDifference"
    }
}
