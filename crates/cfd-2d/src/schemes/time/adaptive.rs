//! Adaptive time stepping for ODE integration
//!
//! This module provides adaptive time step control based on:
//! - CFL condition monitoring (stability-based adaptation)
//! - Richardson extrapolation (accuracy-based adaptation)
//! - Stability monitoring and divergence recovery
//!
//! ## References
//!
//! - Hairer & Wanner (1996): "Solving Ordinary Differential Equations II"
//! - Press et al. (1992): "Numerical Recipes in C"
//! - CFD literature on adaptive time stepping

use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
use std::fmt;

/// Adaptation strategy for time step control
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AdaptationStrategy {
    /// CFL-based adaptation: dt = CFL_target * min(dx/|u|, dy/|v|)
    CFLBased {
        /// Target CFL number (typically 0.5-0.8 for stability)
        cfl_target: f64,
        /// Safety factor for CFL adaptation (typically 0.8-0.9)
        safety_factor: f64,
    },
    /// Error-based adaptation using Richardson extrapolation
    ErrorBased {
        /// Target local truncation error tolerance
        error_tolerance: f64,
        /// Safety factor for error-based adaptation
        safety_factor: f64,
        /// Minimum allowed time step
        dt_min: f64,
        /// Maximum allowed time step
        dt_max: f64,
    },
    /// Combined CFL and error-based adaptation
    Combined {
        /// Target CFL number
        cfl_target: f64,
        /// Error tolerance for accuracy control
        error_tolerance: f64,
        /// Safety factor
        safety_factor: f64,
        /// Minimum time step
        dt_min: f64,
        /// Maximum time step
        dt_max: f64,
    },
}

impl Default for AdaptationStrategy {
    fn default() -> Self {
        Self::CFLBased {
            cfl_target: 0.7,
            safety_factor: 0.8,
        }
    }
}

/// Adaptive time step controller
#[derive(Debug, Clone)]
pub struct AdaptiveController<T: RealField + Copy> {
    /// Current time step
    pub dt_current: T,
    /// Minimum allowed time step
    pub dt_min: T,
    /// Maximum allowed time step
    pub dt_max: T,
    /// Adaptation strategy
    pub strategy: AdaptationStrategy,
    /// Step acceptance counter
    pub steps_accepted: usize,
    /// Step rejection counter
    pub steps_rejected: usize,
    /// Maximum allowed rejections before error
    pub max_rejections: usize,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> AdaptiveController<T> {
    /// Create new adaptive controller with default parameters
    pub fn new(dt_initial: T, strategy: AdaptationStrategy) -> Self {
        let dt_min = T::from_f64(1e-12).unwrap_or_else(T::zero);
        let dt_max = T::from_f64(1e6).unwrap_or_else(T::one);

        Self {
            dt_current: dt_initial,
            dt_min,
            dt_max,
            strategy,
            steps_accepted: 0,
            steps_rejected: 0,
            max_rejections: 10,
        }
    }

    /// Calculate optimal time step based on CFL condition
    ///
    /// # Arguments
    /// * `u_max` - Maximum velocity in x-direction
    /// * `v_max` - Maximum velocity in y-direction
    /// * `dx` - Grid spacing in x-direction
    /// * `dy` - Grid spacing in y-direction
    ///
    /// # Returns
    /// Optimal time step based on CFL condition
    pub fn cfl_based_dt(&self, u_max: T, v_max: T, dx: T, dy: T) -> T {
        let (cfl_target, safety_factor) = match self.strategy {
            AdaptationStrategy::CFLBased {
                cfl_target,
                safety_factor,
            } => (cfl_target, safety_factor),
            AdaptationStrategy::Combined {
                cfl_target,
                safety_factor,
                ..
            } => (cfl_target, safety_factor),
            _ => (0.7, 0.8), // Default values
        };

        // CFL condition: dt ≤ CFL_target * min(dx/|u_max|, dy/|v_max|)
        let dt_cfl_x = if u_max.abs() > T::zero() {
            dx * T::from_f64(cfl_target).unwrap_or_else(T::one) / u_max.abs()
        } else {
            T::from_f64(1e10).unwrap_or_else(T::one) // Large value for zero velocity
        };

        let dt_cfl_y = if v_max.abs() > T::zero() {
            dy * T::from_f64(cfl_target).unwrap_or_else(T::one) / v_max.abs()
        } else {
            T::from_f64(1e10).unwrap_or_else(T::one)
        };

        let dt_cfl = dt_cfl_x.min(dt_cfl_y);

        // Apply safety factor and clamp to bounds
        let dt_adapted = dt_cfl * T::from_f64(safety_factor).unwrap_or_else(T::one);
        dt_adapted.max(self.dt_min).min(self.dt_max)
    }

    /// Estimate local truncation error using Richardson extrapolation
    ///
    /// # Arguments
    /// * `y1` - Solution with step size h
    /// * `y2` - Solution with step size h/2
    /// * `p` - Order of the method
    ///
    /// # Returns
    /// Estimated local truncation error
    pub fn richardson_error(&self, y1: &DVector<T>, y2: &DVector<T>, p: usize) -> T {
        let n = y1.len();
        let mut error_max = T::zero();

        for i in 0..n {
            let diff = (y1[i] - y2[i]).abs();
            let error =
                diff / (T::from_f64(2.0_f64.powi(p as i32)).unwrap_or_else(T::one) - T::one());
            error_max = error_max.max(error);
        }

        error_max
    }

    /// Adapt time step based on error estimate
    ///
    /// # Arguments
    /// * `error_estimate` - Current error estimate
    ///
    /// # Returns
    /// (new_dt, step_accepted)
    pub fn adapt_step(&mut self, error_estimate: T) -> (T, bool) {
        let (error_tolerance, safety_factor) = match self.strategy {
            AdaptationStrategy::ErrorBased {
                error_tolerance,
                safety_factor,
                ..
            } => (error_tolerance, safety_factor),
            AdaptationStrategy::Combined {
                error_tolerance,
                safety_factor,
                ..
            } => (error_tolerance, safety_factor),
            _ => return (self.dt_current, true), // No error-based adaptation
        };

        let error_tolerance_t = T::from_f64(error_tolerance).unwrap_or_else(T::one);

        if error_estimate <= error_tolerance_t {
            // Step accepted - try to increase time step
            self.steps_accepted += 1;

            let factor = T::from_f64(safety_factor).unwrap_or_else(T::one)
                * (error_tolerance_t / error_estimate)
                    .powf(T::from_f64(1.0 / 4.0).unwrap_or_else(T::one));

            let new_dt = (self.dt_current * factor).min(self.dt_max);
            self.dt_current = new_dt;
            (new_dt, true)
        } else {
            // Step rejected - reduce time step
            self.steps_rejected += 1;

            assert!(
                self.steps_rejected <= self.max_rejections,
                "Too many step rejections in adaptive time stepping"
            );

            let factor = T::from_f64(safety_factor).unwrap_or_else(T::one)
                * (error_tolerance_t / error_estimate)
                    .powf(T::from_f64(1.0 / 3.0).unwrap_or_else(T::one));

            let new_dt = (self.dt_current * factor).max(self.dt_min);
            self.dt_current = new_dt;
            (new_dt, false)
        }
    }

    /// Check for numerical divergence
    ///
    /// # Arguments
    /// * `residual` - Current residual norm
    /// * `_residual_prev` - Previous residual norm
    /// * `tolerance` - Divergence tolerance (e.g., 1e10)
    ///
    /// # Returns
    /// true if divergence detected
    pub fn detect_divergence(&self, residual: T, _residual_prev: T, tolerance: T) -> bool {
        // Check if residual is growing rapidly
        if _residual_prev > T::zero() {
            residual / _residual_prev > tolerance
        } else {
            residual > tolerance
        }
    }

    /// Reset controller statistics
    pub fn reset_statistics(&mut self) {
        self.steps_accepted = 0;
        self.steps_rejected = 0;
    }
}

/// Adaptive time integrator combining base integrator with adaptive control
pub struct AdaptiveTimeIntegrator<T: RealField + Copy> {
    /// Base time integrator
    base_integrator: super::TimeIntegrator<T>,
    /// Adaptive controller
    controller: AdaptiveController<T>,
    /// Previous solution for Richardson extrapolation
    y_prev: Option<DVector<T>>,
    /// Previous residual for divergence detection
    _residual_prev: Option<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive> AdaptiveTimeIntegrator<T> {
    /// Create new adaptive time integrator
    pub fn new(scheme: super::TimeScheme, controller: AdaptiveController<T>) -> Self {
        let base_integrator = super::TimeIntegrator::new(scheme);

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
        y: &DVector<T>,
        t: T,
        u_max: T,
        v_max: T,
        dx: T,
        dy: T,
    ) -> (DVector<T>, T, T, bool)
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        // Calculate optimal time step based on CFL condition
        let dt_optimal = self.controller.cfl_based_dt(u_max, v_max, dx, dy);

        // Take step with optimal time step
        let y_new = self.base_integrator.step(&f, y, t, dt_optimal);

        // For CFL-based adaptation, always accept the step
        // (CFL condition ensures stability)
        self.controller.steps_accepted += 1;

        // Update previous solution
        self.y_prev = Some(y.clone());

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
    pub fn step_error_adaptive<F>(&mut self, f: F, y: &DVector<T>, t: T) -> (DVector<T>, T, T, bool)
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let mut attempts = 0;
        let max_attempts = 5;

        loop {
            attempts += 1;
            assert!(
                attempts <= max_attempts,
                "Failed to converge in adaptive time stepping after {max_attempts} attempts"
            );

            let dt_current = self.controller.dt_current;

            // Take full step
            let y_full = self.base_integrator.step(&f, y, t, dt_current);

            // Take two half steps for error estimation
            let y_half = self.base_integrator.step(
                &f,
                y,
                t,
                dt_current / T::from_f64(2.0).unwrap_or_else(T::one),
            );
            let t_half = t + dt_current / T::from_f64(2.0).unwrap_or_else(T::one);
            let y_full_from_half = self.base_integrator.step(
                &f,
                &y_half,
                t_half,
                dt_current / T::from_f64(2.0).unwrap_or_else(T::one),
            );

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
                self.y_prev = Some(y.clone());
                return (y_full, t + dt_current, new_dt, true);
            }
            // Step rejected - try again with smaller step
            continue;
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
        y: &DVector<T>,
        t: T,
        u_max: T,
        v_max: T,
        dx: T,
        dy: T,
    ) -> (DVector<T>, T, T, bool)
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        // First, get CFL-based time step suggestion
        let dt_cfl = self.controller.cfl_based_dt(u_max, v_max, dx, dy);

        // Use minimum of CFL and current adaptive step
        let dt_current = self.controller.dt_current.min(dt_cfl);

        // Now do error-based adaptation with this constrained step
        let mut attempts = 0;
        let max_attempts = 5;

        loop {
            attempts += 1;
            if attempts > max_attempts {
                // Fall back to CFL-only adaptation if error control fails
                return self.step_cfl_adaptive(f, y, t, u_max, v_max, dx, dy);
            }

            // Take full step
            let y_full = self.base_integrator.step(&f, y, t, dt_current);

            // Take two half steps for error estimation
            let y_half = self.base_integrator.step(
                &f,
                y,
                t,
                dt_current / T::from_f64(2.0).unwrap_or_else(T::one),
            );
            let t_half = t + dt_current / T::from_f64(2.0).unwrap_or_else(T::one);
            let y_full_from_half = self.base_integrator.step(
                &f,
                &y_half,
                t_half,
                dt_current / T::from_f64(2.0).unwrap_or_else(T::one),
            );

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
                } => T::from_f64(error_tolerance).unwrap_or_else(T::one),
                _ => T::from_f64(1e-6).unwrap_or_else(T::one), // Default tolerance
            };

            if error_estimate <= error_tolerance {
                // Step accepted
                self.controller.steps_accepted += 1;

                // Try to increase step size for next step
                let factor = T::from_f64(0.9).unwrap_or_else(T::one)
                    * (error_tolerance / error_estimate)
                        .powf(T::from_f64(1.0 / 4.0).unwrap_or_else(T::one));
                self.controller.dt_current = (dt_current * factor).min(self.controller.dt_max);

                self.y_prev = Some(y.clone());
                return (y_full, t + dt_current, self.controller.dt_current, true);
            } else {
                // Step rejected - reduce step size
                self.controller.steps_rejected += 1;
                let factor = T::from_f64(0.9).unwrap_or_else(T::one)
                    * (error_tolerance / error_estimate)
                        .powf(T::from_f64(1.0 / 3.0).unwrap_or_else(T::one));
                self.controller.dt_current = (dt_current * factor).max(self.controller.dt_min);

                // Retry with smaller step
                continue;
            }
        }
    }

    /// Get current time step
    #[must_use]
    pub fn current_dt(&self) -> T {
        self.controller.dt_current
    }

    /// Set current time step
    pub fn set_current_dt(&mut self, dt: T) {
        self.controller.dt_current = dt.max(self.controller.dt_min).min(self.controller.dt_max);
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

impl<T: RealField + Copy> fmt::Display for AdaptiveController<T>
where
    T: fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "AdaptiveController {{ dt: {:.6}, accepted: {}, rejected: {} }}",
            self.dt_current, self.steps_accepted, self.steps_rejected
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::schemes::TimeScheme;
    use approx::assert_relative_eq;

    /// Test ODE: dy/dt = -2y, exact solution: y(t) = y0 * exp(-2t)
    fn test_ode(t: f64, y: &DVector<f64>) -> DVector<f64> {
        let mut f = DVector::zeros(y.len());
        for i in 0..y.len() {
            f[i] = -2.0 * y[i];
        }
        f
    }

    fn exact_solution(t: f64, y0: f64) -> f64 {
        y0 * (-2.0 * t).exp()
    }

    #[test]
    fn test_cfl_based_adaptation() {
        let strategy = AdaptationStrategy::CFLBased {
            cfl_target: 0.7,
            safety_factor: 0.8,
        };

        let mut controller = AdaptiveController::new(0.01, strategy);

        // Test CFL-based dt calculation
        let dt = controller.cfl_based_dt(1.0, 0.5, 0.01, 0.01);
        let expected_dt = 0.01 * 0.7 * 0.8 / 1.0; // CFL limited by u
        assert_relative_eq!(dt, expected_dt, epsilon = 1e-10);

        // Test with zero velocity (should give maximum dt)
        let dt_zero = controller.cfl_based_dt(0.0, 0.0, 0.01, 0.01);
        assert!(dt_zero >= controller.dt_max);
    }

    #[test]
    fn test_error_based_adaptation() {
        let strategy = AdaptationStrategy::ErrorBased {
            error_tolerance: 1e-6,
            safety_factor: 0.9,
            dt_min: 1e-8,
            dt_max: 1.0,
        };

        let mut controller = AdaptiveController::new(0.01, strategy);

        // Test step acceptance (small error)
        let (dt1, accepted1) = controller.adapt_step(1e-8);
        assert!(accepted1);
        assert!(dt1 > 0.01); // Should increase step size

        // Test step rejection (large error)
        let (dt2, accepted2) = controller.adapt_step(1e-3);
        assert!(!accepted2);
        assert!(dt2 < 0.01); // Should decrease step size
    }

    #[test]
    fn test_adaptive_integration_cfl() {
        let strategy = AdaptationStrategy::CFLBased {
            cfl_target: 0.5,
            safety_factor: 0.9,
        };

        let controller = AdaptiveController::new(0.01, strategy);
        let mut integrator = AdaptiveTimeIntegrator::new(TimeScheme::RungeKutta4, controller);

        let y0 = DVector::from_vec(vec![1.0]);
        let t0 = 0.0;

        // Single step with CFL adaptation
        let (y1, t1, dt1, accepted) =
            integrator.step_cfl_adaptive(test_ode, &y0, t0, 1.0, 0.5, 0.01, 0.01);

        assert!(accepted);
        assert!(t1 > t0);
        assert!(dt1 > 0.0);

        // Check that solution is reasonable
        let exact = exact_solution(t1, 1.0);
        assert!((y1[0] - exact).abs() < 0.1); // RK4 should be accurate
    }

    #[test]
    fn test_adaptive_integration_error() {
        let strategy = AdaptationStrategy::ErrorBased {
            error_tolerance: 1e-6,
            safety_factor: 0.9,
            dt_min: 1e-8,
            dt_max: 1.0,
        };

        let controller = AdaptiveController::new(0.01, strategy);
        let mut integrator = AdaptiveTimeIntegrator::new(TimeScheme::RungeKutta4, controller);

        let y0 = DVector::from_vec(vec![1.0]);
        let t0 = 0.0;

        // Single step with error adaptation
        let (y1, t1, dt1, accepted) = integrator.step_error_adaptive(test_ode, &y0, t0);

        assert!(accepted);
        assert!(t1 > t0);
        assert!(dt1 > 0.0);

        // Check that solution is accurate
        let exact = exact_solution(t1, 1.0);
        assert!((y1[0] - exact).abs() < 1e-5); // Should be very accurate
    }
}
