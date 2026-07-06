//! Adaptive time stepping with error control.
//!
//! This module provides adaptive time step controllers that automatically
//! adjust step sizes based on local error estimates to maintain solution
//! accuracy while optimizing computational efficiency.

use super::traits::{
    from_f64, one, state_len, state_norm, state_zeros, zero, EmbeddedMethod, TimeState,
    TimeStepController, TimeStepper,
};
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};

/// Adaptive time stepper using embedded Runge-Kutta methods
pub struct AdaptiveTimeStepper<T: RealField + Copy, M: EmbeddedMethod<T>> {
    /// Base integration method
    method: M,
    /// Time step controller
    controller: StandardController<T>,
    /// Safety factor for step size adjustment
    safety_factor: T,
    /// Minimum allowed time step
    dt_min: T,
    /// Maximum allowed time step
    dt_max: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FloatElement, M: EmbeddedMethod<T>> AdaptiveTimeStepper<T, M> {
    /// Create new adaptive time stepper
    pub fn new(method: M) -> Self {
        Self {
            method,
            controller: StandardController::new(),
            safety_factor: from_f64(0.9),
            dt_min: from_f64(1e-12),
            dt_max: from_f64(1.0),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create with custom parameters
    pub fn with_params(method: M, safety_factor: T, dt_min: T, dt_max: T) -> Self {
        Self {
            method,
            controller: StandardController::new(),
            safety_factor,
            dt_min,
            dt_max,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Take an adaptive step
    pub fn adaptive_step<F>(
        &self,
        f: F,
        t: T,
        u: &TimeState<T>,
        dt_current: T,
        tolerance: T,
    ) -> Result<(TimeState<T>, T, bool)>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        // Take step with error estimation
        let (u_new, error) = self.method.embedded_step(f, t, u, dt_current)?;

        // Compute error norm
        let error_norm = state_norm(&error);

        // Check if step should be accepted
        let accepted = self.controller.accept_step(error_norm, tolerance);

        // Compute optimal new step size
        let dt_new = if accepted {
            let dt_optimal = self
                .controller
                .adapt_step(error_norm, dt_current, tolerance);
            let dt_safe = dt_optimal * self.safety_factor;

            // Clamp to allowed range
            dt_safe.max_scalar(self.dt_min).min_scalar(self.dt_max)
        } else {
            // Reject step, reduce time step
            dt_current * from_f64(0.5)
        };

        Ok((u_new, dt_new, accepted))
    }
}

/// Standard PI controller for time step adaptation
///
/// Based on Hairer & Nørsett's PI controller for ODE solvers
pub struct StandardController<T: RealField + Copy> {
    /// Controller gain for proportional term
    kp: T,
    /// Controller gain for integral term
    ki: T,
}

impl<T: RealField + Copy + FloatElement> Default for StandardController<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FloatElement> StandardController<T> {
    /// Create a new PI controller for adaptive time stepping
    pub fn new() -> Self {
        Self {
            kp: from_f64(0.075),
            ki: from_f64(0.175),
        }
    }
}

impl<T: RealField + Copy + FloatElement> TimeStepController<T> for StandardController<T> {
    fn adapt_step(&self, error_estimate: T, current_dt: T, tolerance: T) -> T {
        if error_estimate <= zero::<T>() {
            return current_dt * from_f64(2.0);
        }

        let error_ratio = tolerance / error_estimate;

        // PI controller: dt_new = dt * error_ratio^(kp + ki * error_ratio)
        // Simplified: dt_new = dt * error_ratio^(0.075 + 0.175 * error_ratio)
        let exponent = self.kp + self.ki * error_ratio;
        let ratio = if exponent > one::<T>() {
            <T as FloatElement>::powf(error_ratio, exponent)
        } else {
            error_ratio
        };

        current_dt * ratio.max_scalar(from_f64(0.1)).min_scalar(from_f64(5.0))
    }
}

/// Embedded Runge-Kutta 4(5) method (Dormand-Prince)
///
/// Provides 4th-order solution and 5th-order error estimate
pub struct DormandPrince54<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for DormandPrince54<T> {
    fn default() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> DormandPrince54<T> {
    /// Create a new Dormand-Prince 5(4) adaptive Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy + FloatElement> TimeStepper<T> for DormandPrince54<T> {
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        let (solution, _error) = self.embedded_step(f, t, u, dt)?;
        Ok(solution)
    }

    fn order(&self) -> usize {
        4
    }

    fn stages(&self) -> usize {
        7
    }

    fn stability_region(&self) -> Option<&str> {
        Some("A-stable for non-stiff problems")
    }
}

impl<T: RealField + Copy + FloatElement> EmbeddedMethod<T> for DormandPrince54<T> {
    fn embedded_step<F>(
        &self,
        f: F,
        t: T,
        u: &TimeState<T>,
        dt: T,
    ) -> Result<(TimeState<T>, TimeState<T>)>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        let n = state_len(u);

        // Dormand-Prince coefficients
        let c = [
            zero::<T>(),
            from_f64(0.2),
            from_f64(0.3),
            from_f64(0.8),
            from_f64(8.0 / 9.0),
            one::<T>(),
            one::<T>(),
        ];

        // Dormand-Prince coefficients (initialize at runtime)
        let mut a = [[zero::<T>(); 7]; 7];
        // Fill lower triangular part of A matrix
        a[1][0] = from_f64(0.2);
        a[2][0] = from_f64(3.0 / 40.0);
        a[2][1] = from_f64(9.0 / 40.0);
        a[3][0] = from_f64(44.0 / 45.0);
        a[3][1] = from_f64(-56.0 / 15.0);
        a[3][2] = from_f64(32.0 / 9.0);
        a[4][0] = from_f64(19372.0 / 6561.0);
        a[4][1] = from_f64(-25360.0 / 2187.0);
        a[4][2] = from_f64(64448.0 / 6561.0);
        a[4][3] = from_f64(-212.0 / 729.0);
        a[5][0] = from_f64(9017.0 / 3168.0);
        a[5][1] = from_f64(-355.0 / 33.0);
        a[5][2] = from_f64(46732.0 / 5247.0);
        a[5][3] = from_f64(49.0 / 176.0);
        a[5][4] = from_f64(-5103.0 / 18656.0);
        a[6][0] = from_f64(35.0 / 384.0);
        a[6][1] = zero::<T>();
        a[6][2] = from_f64(500.0 / 1113.0);
        a[6][3] = from_f64(125.0 / 192.0);
        a[6][4] = from_f64(-2187.0 / 6784.0);
        a[6][5] = from_f64(11.0 / 84.0);

        // 4th-order solution coefficients
        let b4 = [
            from_f64(35.0 / 384.0),
            zero::<T>(),
            from_f64(500.0 / 1113.0),
            from_f64(125.0 / 192.0),
            from_f64(-2187.0 / 6784.0),
            from_f64(11.0 / 84.0),
            zero::<T>(),
        ];

        // 5th-order solution coefficients (for error estimation)
        let b5 = [
            from_f64(5179.0 / 57600.0),
            zero::<T>(),
            from_f64(7571.0 / 16695.0),
            from_f64(393.0 / 640.0),
            from_f64(-92097.0 / 339200.0),
            from_f64(187.0 / 2100.0),
            from_f64(1.0 / 40.0),
        ];

        // Compute stages
        let mut k = vec![state_zeros(n); 7];
        k[0] = f(t, u)?;

        for stage in 1..7 {
            let mut u_stage = state_zeros(n);
            for i in 0..n {
                u_stage[i] = u[i];
                for j in 0..stage {
                    u_stage[i] += dt * a[stage][j] * k[j][i];
                }
            }
            let t_stage = t + c[stage] * dt;
            k[stage] = f(t_stage, &u_stage)?;
        }

        // Compute 4th-order solution
        let mut u4 = state_zeros(n);
        for i in 0..n {
            for stage in 0..7 {
                u4[i] += b4[stage] * k[stage][i];
            }
            u4[i] = u[i] + dt * u4[i];
        }

        // Compute 5th-order solution for error estimation
        let mut u5 = state_zeros(n);
        for i in 0..n {
            for stage in 0..7 {
                u5[i] += b5[stage] * k[stage][i];
            }
            u5[i] = u[i] + dt * u5[i];
        }

        // Error estimate: difference between 4th and 5th order solutions
        let mut error = state_zeros(n);
        for i in 0..n {
            error[i] = u5[i] - u4[i];
        }

        Ok((u4, error))
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::time_stepping::traits::{state_from_vec, state_neg};
    use approx::assert_relative_eq;

    // Test function: du/dt = -u (exponential decay)
    fn exponential_decay(_t: f64, u: &TimeState<f64>) -> Result<TimeState<f64>> {
        Ok(state_neg(u))
    }

    #[test]
    fn test_dormand_prince_properties() {
        let dp = DormandPrince54::<f64>::new();
        assert_eq!(dp.order(), 4);
        assert_eq!(dp.stages(), 7);
    }

    #[test]
    fn test_adaptive_stepper() {
        let method = DormandPrince54::new();
        let stepper = AdaptiveTimeStepper::new(method);

        let u0 = state_from_vec(vec![1.0]);
        let dt_initial = 0.1;
        let tolerance = 1e-6;
        let t = 0.0;

        let (u_new, dt_new, accepted) = stepper
            .adaptive_step(exponential_decay, t, &u0, dt_initial, tolerance)
            .unwrap();

        assert!(accepted);
        assert!(dt_new > 0.0);

        // Should produce reasonable result
        let u_analytical = 1.0 * (-dt_initial).exp();
        assert_relative_eq!(u_new[0], u_analytical, epsilon = 1e-4);
    }
}
