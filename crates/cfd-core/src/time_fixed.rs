//! Time integration module for CFD simulations
//!
//! Provides various time integration schemes including explicit and implicit methods.

use crate::error::{Error, Result};
use nalgebra::{DVector, RealField};
use serde::{Deserialize, Serialize};

/// Time integration scheme trait
pub trait TimeIntegrator<T: RealField + Copy>: Send + Sync {
    /// State type for the integrator
    type State;

    /// Perform a single time step
    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Get the integrator name
    fn name(&self) -> &str;
}

/// Forward Euler (explicit) time integration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForwardEuler;

impl<T: RealField + Copy> TimeIntegrator<T> for ForwardEuler {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let f_val = f(t, state);
        *state += f_val * dt;
        Ok(())
    }

    fn order(&self) -> usize {
        1
    }

    fn name(&self) -> &str {
        "Forward Euler"
    }
}

/// Backward Euler (implicit) time integration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BackwardEuler<T: RealField + Copy> {
    /// Maximum iterations for fixed-point iteration
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
}

impl<T: RealField + Copy> Default for BackwardEuler<T> {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::default_epsilon()),
        }
    }
}

impl<T: RealField + Copy> TimeIntegrator<T> for BackwardEuler<T> {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Implement implicit Backward Euler using fixed-point iteration
        // Solve: y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})

        let y_old = state.clone();
        let mut y_new = state.clone();
        let t_new = t + dt;

        // Check for valid iteration count
        if self.max_iterations == 0 {
            return Err(Error::InvalidConfiguration(
                "BackwardEuler requires max_iterations > 0".to_string()
            ));
        }

        // Fixed-point iteration: y_{n+1}^{k+1} = y_n + dt * f(t_{n+1}, y_{n+1}^k)
        for iteration in 0..self.max_iterations {
            let f_val = f(t_new, &y_new);
            let y_next = &y_old + f_val * dt;

            // Check convergence
            let error = (&y_next - &y_new).norm();
            if error < self.tolerance {
                *state = y_next;
                return Ok(());
            }

            y_new = y_next;
        }

        Err(Error::ConvergenceError(format!(
            "Backward Euler did not converge after {} iterations",
            self.max_iterations
        )))
    }

    fn order(&self) -> usize {
        1
    }

    fn name(&self) -> &str {
        "Backward Euler"
    }
}

/// Crank-Nicolson (implicit) time integration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CrankNicolson<T: RealField + Copy> {
    /// Maximum iterations for fixed-point iteration
    pub max_iterations: usize,
    /// Convergence tolerance
    pub tolerance: T,
}

impl<T: RealField + Copy> Default for CrankNicolson<T> {
    fn default() -> Self {
        Self {
            max_iterations: 100,
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::default_epsilon()),
        }
    }
}

impl<T: RealField + Copy> TimeIntegrator<T> for CrankNicolson<T> {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Crank-Nicolson: y_{n+1} = y_n + (dt/2) * [f(t_n, y_n) + f(t_{n+1}, y_{n+1})]
        
        let y_old = state.clone();
        let mut y_new = state.clone();
        let t_new = t + dt;
        let half_dt = dt / (T::one() + T::one());
        
        // Evaluate f at current time
        let f_old = f(t, &y_old);
        
        // Fixed-point iteration
        for _ in 0..self.max_iterations {
            let f_new = f(t_new, &y_new);
            let y_next = &y_old + (&f_old + &f_new) * half_dt;
            
            // Check convergence
            let error = (&y_next - &y_new).norm();
            if error < self.tolerance {
                *state = y_next;
                return Ok(());
            }
            
            y_new = y_next;
        }
        
        Err(Error::ConvergenceError(format!(
            "Crank-Nicolson did not converge after {} iterations",
            self.max_iterations
        )))
    }

    fn order(&self) -> usize {
        2
    }

    fn name(&self) -> &str {
        "Crank-Nicolson"
    }
}

/// Runge-Kutta 4th order time integration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RungeKutta4;

impl<T: RealField + Copy> TimeIntegrator<T> for RungeKutta4 {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let two = T::one() + T::one();
        let six = two * (two + T::one());
        
        // RK4 stages
        let k1 = f(t, state);
        let y_temp = state.clone() + &k1 * (dt / two);
        
        let k2 = f(t + dt / two, &y_temp);
        let y_temp = state.clone() + &k2 * (dt / two);
        
        let k3 = f(t + dt / two, &y_temp);
        let y_temp = state.clone() + &k3 * dt;
        
        let k4 = f(t + dt, &y_temp);
        
        // Update state
        *state += (k1 + k2 * two + k3 * two + k4) * (dt / six);
        
        Ok(())
    }

    fn order(&self) -> usize {
        4
    }

    fn name(&self) -> &str {
        "Runge-Kutta 4"
    }
}

/// Adaptive time stepping controller
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AdaptiveTimeStepController<T: RealField + Copy> {
    /// Minimum allowed time step
    pub dt_min: T,
    /// Maximum allowed time step
    pub dt_max: T,
    /// Target relative error
    pub target_error: T,
    /// Safety factor for time step adjustment
    pub safety_factor: T,
    /// Maximum factor for time step increase
    pub max_increase_factor: T,
    /// Maximum factor for time step decrease
    pub max_decrease_factor: T,
}

impl<T: RealField + Copy> Default for AdaptiveTimeStepController<T> {
    fn default() -> Self {
        Self {
            dt_min: T::from_f64(1e-10).unwrap_or_else(|| T::default_epsilon()),
            dt_max: T::from_f64(1.0).unwrap_or_else(T::one),
            target_error: T::from_f64(1e-6).unwrap_or_else(|| T::default_epsilon() * T::from_f64(100.0).unwrap()),
            safety_factor: T::from_f64(0.9).unwrap_or_else(|| T::one() - T::from_f64(0.1).unwrap()),
            max_increase_factor: T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one()),
            max_decrease_factor: T::from_f64(0.1).unwrap_or_else(|| T::from_f64(0.1).unwrap()),
        }
    }
}

impl<T: RealField + Copy> AdaptiveTimeStepController<T> {
    /// Compute new time step based on error estimate
    pub fn compute_new_dt(&self, current_dt: T, error: T, order: usize) -> T {
        if error < T::default_epsilon() {
            // Error is essentially zero, increase time step
            return (current_dt * self.max_increase_factor).min(self.dt_max);
        }
        
        // Compute optimal time step using error estimate
        let order_t = T::from_usize(order).unwrap_or_else(T::one);
        let factor = self.safety_factor * (self.target_error / error).powf(T::one() / (order_t + T::one()));
        
        // Limit the change factor
        let factor = factor.max(self.max_decrease_factor).min(self.max_increase_factor);
        
        // Apply limits
        (current_dt * factor).max(self.dt_min).min(self.dt_max)
    }
    
    /// Check if time step should be rejected
    pub fn should_reject(&self, error: T) -> bool {
        error > self.target_error
    }
}

/// Time integration configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeIntegrationConfig<T: RealField + Copy> {
    /// Integration scheme
    pub scheme: TimeScheme,
    /// Time step size
    pub dt: T,
    /// Start time
    pub t_start: T,
    /// End time
    pub t_end: T,
    /// Enable adaptive time stepping
    pub adaptive: bool,
    /// Adaptive controller configuration
    pub adaptive_controller: Option<AdaptiveTimeStepController<T>>,
}

/// Available time integration schemes
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TimeScheme {
    /// Forward Euler
    ForwardEuler,
    /// Backward Euler
    BackwardEuler,
    /// Crank-Nicolson
    CrankNicolson,
    /// Runge-Kutta 4th order
    RungeKutta4,
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;
    
    #[test]
    fn test_forward_euler() {
        let integrator = ForwardEuler;
        let mut state = DVector::from_vec(vec![1.0, 0.0]);
        
        // Simple harmonic oscillator: dy/dt = [v, -y]
        let f = |_t: f64, y: &DVector<f64>| -> DVector<f64> {
            DVector::from_vec(vec![y[1], -y[0]])
        };
        
        let dt = 0.01;
        let result = integrator.step(&mut state, 0.0, dt, f);
        
        assert!(result.is_ok());
        assert!((state[0] - 1.0).abs() < 0.01);
        assert!((state[1] + 0.01).abs() < 0.01);
    }
    
    #[test]
    fn test_runge_kutta4() {
        let integrator = RungeKutta4;
        let mut state = DVector::from_vec(vec![1.0, 0.0]);
        
        // Simple harmonic oscillator
        let f = |_t: f64, y: &DVector<f64>| -> DVector<f64> {
            DVector::from_vec(vec![y[1], -y[0]])
        };
        
        let dt = 0.1;
        let result = integrator.step(&mut state, 0.0, dt, f);
        
        assert!(result.is_ok());
        // RK4 should be much more accurate than Forward Euler
        assert!((state[0] - 0.995).abs() < 0.001);
    }
    
    #[test]
    fn test_adaptive_controller() {
        let controller = AdaptiveTimeStepController::<f64>::default();
        
        // Test time step increase for small error
        let new_dt = controller.compute_new_dt(0.01, 1e-8, 4);
        assert!(new_dt > 0.01);
        
        // Test time step decrease for large error
        let new_dt = controller.compute_new_dt(0.01, 1e-3, 4);
        assert!(new_dt < 0.01);
        
        // Test rejection
        assert!(controller.should_reject(1e-3));
        assert!(!controller.should_reject(1e-8));
    }
}