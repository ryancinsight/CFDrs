//! Adaptive time integrator combining base scheme with adaptive step control
//!
//! Wraps a [`super::super::TimeIntegrator`] with an [`AdaptiveController`] to
//! provide CFL-based, error-based, or combined adaptive time stepping.

use super::{AdaptationStrategy, AdaptiveController};
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use eunomia::{FloatElement, NumericElement};

use super::super::StateVector;

/// Adaptive time integrator combining base integrator with adaptive control
pub struct AdaptiveTimeIntegrator<T: Cfd2dScalar + Copy> {
    /// Base time integrator
    base_integrator: super::super::TimeIntegrator<T>,
    /// Adaptive controller
    controller: AdaptiveController<T>,
    /// Previous solution for Richardson extrapolation
    y_prev: Option<StateVector<T>>,
    /// Previous residual for divergence detection
    _residual_prev: Option<T>,
}

impl<T: Cfd2dScalar + Copy + FloatElement> AdaptiveTimeIntegrator<T> {
    /// Update the y_prev buffer used by Richardson-extrapolation error estimators.
    ///
    /// AUDIT: allocates ONCE on the first call (when the `Option` is `None`),
    /// then reuses the existing `Array1` buffer on every subsequent call via
    /// in-place memcpy (`prev.assign(y)`). Saves 1 Vec-clone allocation per
    /// call after the first. Functionally identical to a `mem::take +
    /// copy_from_slice` pattern, just with a tighter API.
    fn store_y_prev(&mut self, y: &StateVector<T>) {
        if let Some(prev) = self.y_prev.as_mut() {
            prev.assign(y);
        } else {
            self.y_prev = Some(y.clone());
        }
    }

    /// Create new adaptive time integrator
    pub fn new(scheme: super::super::TimeScheme, controller: AdaptiveController<T>) -> Self {
        let base_integrator = super::super::TimeIntegrator::new(scheme);

        Self {
            base_integrator,
            controller,
            y_prev: None,
            _residual_prev: None,
        }
    }

    /// Perform adaptive time step with CFL-based adaptation
    ///
    /// # Arguments
    /// * `f` - Right-hand side function
    /// * `y` - Current solution
    /// * `t` - Current time
    /// * `u_max` - Maximum velocity in x-direction
    /// * `v_max` - Maximum velocity in y-direction
    /// * `dx` - Grid spacing in x-direction
    /// * `dy` - Grid spacing in y-direction
    ///
    /// # Returns
    /// (solution, new_time, new_dt, step_accepted)
    pub fn step_cfl_adaptive<F>(
        &mut self,
        f: F,
        y: &StateVector<T>,
        t: T,
        u_max: T,
        v_max: T,
        dx: T,
        dy: T,
    ) -> (StateVector<T>, T, T, bool)
    where
        F: Fn(T, &StateVector<T>) -> StateVector<T>,
    {
        // Calculate optimal time step based on CFL condition
        let dt_optimal = self.controller.cfl_based_dt(u_max, v_max, dx, dy);

        // Take step with optimal time step
        let y_new = self.base_integrator.step(&f, y, t, dt_optimal);

        // For CFL-based adaptation, always accept the step
        // (CFL condition ensures stability)
        self.controller.steps_accepted += 1;
        self.store_y_prev(y);

        (y_new, t + dt_optimal, dt_optimal, true)
    }

    /// Perform adaptive time step with error-based adaptation
    ///
    /// # Arguments
    /// * `f` - Right-hand side function
    /// * `y` - Current solution
    /// * `t` - Current time
    ///
    /// # Returns
    /// (solution, new_time, new_dt, step_accepted)
    pub fn step_error_adaptive<F>(
        &mut self,
        f: F,
        y: &StateVector<T>,
        t: T,
    ) -> (StateVector<T>, T, T, bool)
    where
        F: Fn(T, &StateVector<T>) -> StateVector<T>,
    {
        let mut attempts = 0;
        let max_attempts = 5;

        loop {
            attempts += 1;
            if attempts > max_attempts {
                let dt_fallback = <T as NumericElement>::max_scalar(
                    self.controller.dt_current,
                    self.controller.dt_min,
                );
                let y_fallback = self.base_integrator.step(&f, y, t, dt_fallback);
                self.controller.dt_current = dt_fallback;
                self.controller.steps_rejected += 1;
                self.controller.steps_accepted += 1;
                self.store_y_prev(y);
                return (y_fallback, t + dt_fallback, dt_fallback, true);
            }

            let dt_current = self.controller.dt_current;

            // Take full step
            let y_full = self.base_integrator.step(&f, y, t, dt_current);

            // Take two half steps for error estimation
            let two = scalar::from_f64::<T>(2.0);
            let half_dt = dt_current / two;
            let y_half = self.base_integrator.step(&f, y, t, half_dt);
            let t_half = t + half_dt;
            let y_full_from_half = self.base_integrator.step(&f, &y_half, t_half, half_dt);

            // Estimate error using Richardson extrapolation
            let error_estimate = self.controller.richardson_error(
                &y_full,
                &y_full_from_half,
                self.base_integrator.order(),
            );

            // Adapt step size based on error
            let (new_dt, accepted) = self.controller.adapt_step(error_estimate);

            if accepted {
                // Step accepted
                self.store_y_prev(y);
                return (y_full, t + dt_current, new_dt, true);
            }
            // Step rejected - try again with smaller step
        }
    }

    /// Perform adaptive time step with combined CFL and error control
    ///
    /// # Arguments
    /// * `f` - Right-hand side function
    /// * `y` - Current solution
    /// * `t` - Current time
    /// * `u_max` - Maximum velocity in x-direction
    /// * `v_max` - Maximum velocity in y-direction
    /// * `dx` - Grid spacing in x-direction
    /// * `dy` - Grid spacing in y-direction
    ///
    /// # Returns
    /// (solution, new_time, new_dt, step_accepted)
    pub fn step_combined_adaptive<F>(
        &mut self,
        f: F,
        y: &StateVector<T>,
        t: T,
        u_max: T,
        v_max: T,
        dx: T,
        dy: T,
    ) -> (StateVector<T>, T, T, bool)
    where
        F: Fn(T, &StateVector<T>) -> StateVector<T>,
    {
        // First, get CFL-based time step suggestion
        let dt_cfl = self.controller.cfl_based_dt(u_max, v_max, dx, dy);

        // Use minimum of CFL and current adaptive step
        let dt_current = <T as NumericElement>::min_scalar(self.controller.dt_current, dt_cfl);

        // Now do error-based adaptation with this constrained step
        let mut attempts = 0;
        let max_attempts = 5;

        loop {
            attempts += 1;
            if attempts > max_attempts {
                return self.step_cfl_adaptive(f, y, t, u_max, v_max, dx, dy);
            }

            // Take full step
            let y_full = self.base_integrator.step(&f, y, t, dt_current);

            // Take two half steps for error estimation
            let two = scalar::from_f64::<T>(2.0);
            let half_dt = dt_current / two;
            let y_half = self.base_integrator.step(&f, y, t, half_dt);
            let t_half = t + half_dt;
            let y_full_from_half = self.base_integrator.step(&f, &y_half, t_half, half_dt);

            // Estimate error using Richardson extrapolation
            let error_estimate = self.controller.richardson_error(
                &y_full,
                &y_full_from_half,
                self.base_integrator.order(),
            );

            // Check if error is acceptable
            let error_tolerance = match self.controller.strategy {
                AdaptationStrategy::Combined {
                    error_tolerance, ..
                } => scalar::from_f64::<T>(error_tolerance),
                _ => scalar::from_f64::<T>(1e-6), // Default tolerance
            };

            if error_estimate <= error_tolerance {
                // Step accepted
                self.controller.steps_accepted += 1;

                // Try to increase step size for next step
                let factor = scalar::from_f64::<T>(0.9)
                    * <T as FloatElement>::powf(
                        error_tolerance / error_estimate,
                        scalar::from_f64::<T>(1.0 / 4.0),
                    );                    self.controller.dt_current =
                        <T as NumericElement>::min_scalar(dt_current * factor, self.controller.dt_max);

                    self.store_y_prev(y);
                    return (y_full, t + dt_current, self.controller.dt_current, true);
                }

            // Step rejected - reduce step size
            self.controller.steps_rejected += 1;
            let factor = scalar::from_f64::<T>(0.9)
                * <T as FloatElement>::powf(
                    error_tolerance / error_estimate,
                    scalar::from_f64::<T>(1.0 / 3.0),
                );
            self.controller.dt_current =
                <T as NumericElement>::max_scalar(dt_current * factor, self.controller.dt_min);

            // Retry with smaller step
        }
    }

    /// Get current time step
    #[must_use]
    pub fn current_dt(&self) -> T {
        self.controller.dt_current
    }

    /// Set current time step
    pub fn set_current_dt(&mut self, dt: T) {
        self.controller.dt_current = <T as NumericElement>::min_scalar(
            <T as NumericElement>::max_scalar(dt, self.controller.dt_min),
            self.controller.dt_max,
        );
    }

    /// Adapt time step based on error estimate
    pub fn adapt_step(&mut self, error_estimate: T) -> (T, bool) {
        self.controller.adapt_step(error_estimate)
    }

    /// Get controller statistics
    pub fn statistics(&self) -> (usize, usize) {
        (
            self.controller.steps_accepted,
            self.controller.steps_rejected,
        )
    }

    /// Reset controller statistics
    pub fn reset_statistics(&mut self) {
        self.controller.reset_statistics();
    }
}
