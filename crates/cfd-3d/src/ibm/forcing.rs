//! Forcing methods for IBM

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Forcing method for IBM
pub trait ForcingMethod<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Calculate forcing term
    fn calculate_force(
        &self,
        desired_velocity: &Vector3<T>,
        current_velocity: &Vector3<T>,
        dt: T,
    ) -> Vector3<T>;

    /// Update forcing based on feedback
    fn update(&mut self, error: &Vector3<T>);
}

/// Direct forcing method implementing the exact Fadlun et al. (2000) formula.
///
/// # Theorem (Direct Forcing, Fadlun et al. 2000, Eq. 8)
///
/// The immersed-boundary body force that enforces no-slip to O(Δt) accuracy is:
///
///   f_ib = (u_desired − u*) / Δt
///
/// There are no free dimensionless parameters. Any scale factor ≠ 1 violates the
/// O(Δt) exactness proof.
pub struct DirectForcing;

impl DirectForcing {
    /// Create a new direct forcing operator (no free parameters)
    #[must_use]
    pub fn new() -> Self {
        Self
    }
}

impl Default for DirectForcing {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> ForcingMethod<T>
    for DirectForcing
{
    fn calculate_force(
        &self,
        desired_velocity: &Vector3<T>,
        current_velocity: &Vector3<T>,
        dt: T,
    ) -> Vector3<T> {
        // Exact direct-forcing: f = (u_desired − u*) / Δt
        if dt > T::zero() {
            (desired_velocity - current_velocity) / dt
        } else {
            Vector3::zeros()
        }
    }

    fn update(&mut self, _error: &Vector3<T>) {
        // Direct forcing has no internal state
    }
}


/// Feedback forcing method with proportional-integral control
pub struct FeedbackForcing<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    kp: T, // Proportional gain
    ki: T, // Integral gain
    integral: Vector3<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> FeedbackForcing<T> {
    /// Create a new feedback forcing method
    pub fn new(kp: T, ki: T) -> Self {
        Self {
            kp,
            ki,
            integral: Vector3::zeros(),
        }
    }

    /// Reset integral term
    pub fn reset(&mut self) {
        self.integral = Vector3::zeros();
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> ForcingMethod<T> for FeedbackForcing<T> {
    fn calculate_force(
        &self,
        desired_velocity: &Vector3<T>,
        current_velocity: &Vector3<T>,
        _dt: T,
    ) -> Vector3<T> {
        let error = desired_velocity - current_velocity;

        // PI control
        error * self.kp + self.integral * self.ki
    }

    fn update(&mut self, error: &Vector3<T>) {
        self.integral += error;
    }
}
