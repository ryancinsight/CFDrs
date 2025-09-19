//! Forcing methods for IBM

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Forcing method for IBM
pub trait ForcingMethod<T: RealField + Copy> {
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

/// Direct forcing method
pub struct DirectForcing<T: RealField + Copy> {
    scale: T,
}

impl<T: RealField + FromPrimitive + Copy> DirectForcing<T> {
    /// Create a new direct forcing method
    pub fn new(scale: T) -> Self {
        Self { scale }
    }
}

impl<T: RealField + FromPrimitive + Copy> ForcingMethod<T> for DirectForcing<T> {
    fn calculate_force(
        &self,
        desired_velocity: &Vector3<T>,
        current_velocity: &Vector3<T>,
        dt: T,
    ) -> Vector3<T> {
        if dt > T::zero() {
            (desired_velocity - current_velocity) * self.scale / dt
        } else {
            Vector3::zeros()
        }
    }

    fn update(&mut self, _error: &Vector3<T>) {
        // Direct forcing doesn't need updates
    }
}

/// Feedback forcing method with proportional-integral control
pub struct FeedbackForcing<T: RealField + Copy> {
    kp: T, // Proportional gain
    ki: T, // Integral gain
    integral: Vector3<T>,
}

impl<T: RealField + FromPrimitive + Copy> FeedbackForcing<T> {
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

impl<T: RealField + FromPrimitive + Copy> ForcingMethod<T> for FeedbackForcing<T> {
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
