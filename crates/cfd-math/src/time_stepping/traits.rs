//! Traits for time-stepping methods.

use cfd_core::error::Result;
use nalgebra::{DVector, RealField};

/// Time-stepping method for ODE integration du/dt = f(t,u)
pub trait TimeStepper<T: RealField + Copy> {
    /// Advance solution by one time step
    ///
    /// # Arguments
    /// * `f` - Right-hand side function f(t, u)
    /// * `t` - Current time
    /// * `u` - Current solution vector
    /// * `dt` - Time step size
    ///
    /// # Returns
    /// New solution vector at t + dt
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>;

    /// Get the order of accuracy of this method
    fn order(&self) -> usize;

    /// Get the number of stages (function evaluations per step)
    fn stages(&self) -> usize;

    /// Check if this is an explicit method
    fn is_explicit(&self) -> bool {
        true // Default to explicit
    }

    /// Get method stability region information
    fn stability_region(&self) -> Option<&str> {
        None
    }
}

/// Controller for adaptive time stepping
pub trait TimeStepController<T: RealField + Copy> {
    /// Compute optimal time step based on error estimate
    ///
    /// # Arguments
    /// * `error_estimate` - Estimated local truncation error
    /// * `current_dt` - Current time step size
    /// * `tolerance` - Desired error tolerance
    ///
    /// # Returns
    /// Recommended new time step size
    fn adapt_step(&self, error_estimate: T, current_dt: T, tolerance: T) -> T;

    /// Check if solution should be accepted
    ///
    /// # Arguments
    /// * `error_estimate` - Estimated local truncation error
    /// * `tolerance` - Desired error tolerance
    ///
    /// # Returns
    /// true if step should be accepted, false if it should be rejected
    fn accept_step(&self, error_estimate: T, tolerance: T) -> bool {
        error_estimate <= tolerance
    }
}

/// Embedded Runge-Kutta method providing error estimation
pub trait EmbeddedMethod<T: RealField + Copy>: TimeStepper<T> {
    /// Take a step and return both solution and error estimate
    ///
    /// # Arguments
    /// * `f` - Right-hand side function f(t, u)
    /// * `t` - Current time
    /// * `u` - Current solution vector
    /// * `dt` - Time step size
    ///
    /// # Returns
    /// Tuple of (solution, error_estimate)
    fn embedded_step<F>(
        &self,
        f: F,
        t: T,
        u: &DVector<T>,
        dt: T,
    ) -> Result<(DVector<T>, DVector<T>)>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>;
}
