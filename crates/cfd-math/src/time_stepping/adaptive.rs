//! Adaptive time stepping with error control.
//!
//! This module provides adaptive time step controllers that automatically
//! adjust step sizes based on local error estimates to maintain solution
//! accuracy while optimizing computational efficiency.

use super::traits::{EmbeddedMethod, TimeStepController, TimeStepper};
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};

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

impl<T: RealField + Copy, M: EmbeddedMethod<T>> AdaptiveTimeStepper<T, M> {
    /// Create new adaptive time stepper
    pub fn new(method: M) -> Self {
        Self {
            method,
            controller: StandardController::new(),
            safety_factor: T::from_f64(0.9).unwrap(),
            dt_min: T::from_f64(1e-12).unwrap(),
            dt_max: T::from_f64(1.0).unwrap(),
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
        u: &DVector<T>,
        dt_current: T,
        tolerance: T,
    ) -> Result<(DVector<T>, T, bool)>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        // Take step with error estimation
        let (u_new, error) = self.method.embedded_step(f, t, u, dt_current)?;

        // Compute error norm
        let error_norm = error.norm();

        // Check if step should be accepted
        let accepted = self.controller.accept_step(error_norm, tolerance);

        // Compute optimal new step size
        let dt_new = if accepted {
            let dt_optimal = self
                .controller
                .adapt_step(error_norm, dt_current, tolerance);
            let dt_safe = dt_optimal * self.safety_factor;

            // Clamp to allowed range
            dt_safe.max(self.dt_min).min(self.dt_max)
        } else {
            // Reject step, reduce time step
            dt_current * T::from_f64(0.5).unwrap()
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
    /// Previous error ratio (for integral term)
    #[allow(dead_code)]
    prev_error_ratio: Option<T>,
}

impl<T: RealField + Copy> Default for StandardController<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> StandardController<T> {
    /// Create a new PI controller for adaptive time stepping
    pub fn new() -> Self {
        Self {
            kp: T::from_f64(0.075).unwrap(),
            ki: T::from_f64(0.175).unwrap(),
            prev_error_ratio: None,
        }
    }
}

impl<T: RealField + Copy> TimeStepController<T> for StandardController<T> {
    fn adapt_step(&self, error_estimate: T, current_dt: T, tolerance: T) -> T {
        if error_estimate <= T::zero() {
            return current_dt * T::from_f64(2.0).unwrap();
        }

        let error_ratio = tolerance / error_estimate;

        // PI controller: dt_new = dt * error_ratio^(kp + ki * error_ratio)
        // Simplified: dt_new = dt * error_ratio^(0.075 + 0.175 * error_ratio)
        let exponent = self.kp + self.ki * error_ratio;
        let ratio = if exponent > T::from_f64(1.0).unwrap() {
            error_ratio.powf(exponent)
        } else {
            error_ratio
        };

        current_dt
            * ratio
                .max(T::from_f64(0.1).unwrap())
                .min(T::from_f64(5.0).unwrap())
    }
}

/// Embedded Runge-Kutta 4(5) method (Dormand-Prince)
///
/// Provides 4th-order solution and 5th-order error estimate
pub struct DormandPrince54<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Default for DormandPrince54<T> {
    fn default() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField> DormandPrince54<T> {
    /// Create a new Dormand-Prince 5(4) adaptive Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy> TimeStepper<T> for DormandPrince54<T> {
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
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

impl<T: RealField + Copy> EmbeddedMethod<T> for DormandPrince54<T> {
    fn embedded_step<F>(
        &self,
        f: F,
        t: T,
        u: &DVector<T>,
        dt: T,
    ) -> Result<(DVector<T>, DVector<T>)>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        let n = u.len();

        // Dormand-Prince coefficients
        let c = [
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.2).unwrap(),
            T::from_f64(0.3).unwrap(),
            T::from_f64(0.8).unwrap(),
            T::from_f64(8.0 / 9.0).unwrap(),
            T::from_f64(1.0).unwrap(),
            T::from_f64(1.0).unwrap(),
        ];

        // Dormand-Prince coefficients (initialize at runtime)
        let mut a = [[T::zero(); 7]; 7];
        // Fill lower triangular part of A matrix
        a[1][0] = T::from_f64(0.2).unwrap();
        a[2][0] = T::from_f64(3.0 / 40.0).unwrap();
        a[2][1] = T::from_f64(9.0 / 40.0).unwrap();
        a[3][0] = T::from_f64(44.0 / 45.0).unwrap();
        a[3][1] = T::from_f64(-56.0 / 15.0).unwrap();
        a[3][2] = T::from_f64(32.0 / 9.0).unwrap();
        a[4][0] = T::from_f64(19372.0 / 6561.0).unwrap();
        a[4][1] = T::from_f64(-25360.0 / 2187.0).unwrap();
        a[4][2] = T::from_f64(64448.0 / 6561.0).unwrap();
        a[4][3] = T::from_f64(-212.0 / 729.0).unwrap();
        a[5][0] = T::from_f64(9017.0 / 3168.0).unwrap();
        a[5][1] = T::from_f64(-355.0 / 33.0).unwrap();
        a[5][2] = T::from_f64(46732.0 / 5247.0).unwrap();
        a[5][3] = T::from_f64(49.0 / 176.0).unwrap();
        a[5][4] = T::from_f64(-5103.0 / 18656.0).unwrap();
        a[6][0] = T::from_f64(35.0 / 384.0).unwrap();
        a[6][1] = T::zero();
        a[6][2] = T::from_f64(500.0 / 1113.0).unwrap();
        a[6][3] = T::from_f64(125.0 / 192.0).unwrap();
        a[6][4] = T::from_f64(-2187.0 / 6784.0).unwrap();
        a[6][5] = T::from_f64(11.0 / 84.0).unwrap();

        // 4th-order solution coefficients
        let b4 = [
            T::from_f64(35.0 / 384.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(500.0 / 1113.0).unwrap(),
            T::from_f64(125.0 / 192.0).unwrap(),
            T::from_f64(-2187.0 / 6784.0).unwrap(),
            T::from_f64(11.0 / 84.0).unwrap(),
            T::from_f64(0.0).unwrap(),
        ];

        // 5th-order solution coefficients (for error estimation)
        let b5 = [
            T::from_f64(5179.0 / 57600.0).unwrap(),
            T::from_f64(0.0).unwrap(),
            T::from_f64(7571.0 / 16695.0).unwrap(),
            T::from_f64(393.0 / 640.0).unwrap(),
            T::from_f64(-92097.0 / 339200.0).unwrap(),
            T::from_f64(187.0 / 2100.0).unwrap(),
            T::from_f64(1.0 / 40.0).unwrap(),
        ];

        // Compute stages
        let mut k = vec![DVector::zeros(n); 7];
        k[0] = f(t, u)?;

        for stage in 1..7 {
            let mut u_stage = DVector::zeros(n);
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
        let mut u4 = DVector::zeros(n);
        for i in 0..n {
            for stage in 0..7 {
                u4[i] += b4[stage] * k[stage][i];
            }
            u4[i] = u[i] + dt * u4[i];
        }

        // Compute 5th-order solution for error estimation
        let mut u5 = DVector::zeros(n);
        for i in 0..n {
            for stage in 0..7 {
                u5[i] += b5[stage] * k[stage][i];
            }
            u5[i] = u[i] + dt * u5[i];
        }

        // Error estimate: difference between 4th and 5th order solutions
        let error = &u5 - &u4;

        Ok((u4, error))
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use approx::assert_relative_eq;

    // Test function: du/dt = -u (exponential decay)
    fn exponential_decay(_t: f64, u: &DVector<f64>) -> Result<DVector<f64>> {
        Ok(-u.clone())
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

        let u0 = DVector::from_vec(vec![1.0]);
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
