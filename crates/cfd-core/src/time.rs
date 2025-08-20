//! Time integration schemes.

use crate::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::cast::FromPrimitive;
use num_traits::Float;


/// Trait for time integration schemes
pub trait TimeIntegrator<T: RealField>: Send + Sync {
    /// State type
    type State;

    /// Perform one time step
    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Check if the scheme is explicit
    fn is_explicit(&self) -> bool;
}

/// Forward Euler (explicit) time integration
pub struct ForwardEuler;

impl<T: RealField> TimeIntegrator<T> for ForwardEuler {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let derivative = f(t, state);
        state.axpy(dt, &derivative, T::one());
        Ok(())
    }

    fn order(&self) -> usize {
        1
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

/// Runge-Kutta 2nd order (RK2) time integration
pub struct RungeKutta2;

impl<T: RealField + FromPrimitive> TimeIntegrator<T> for RungeKutta2 {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation)
        })?;

        let k1 = f(t.clone(), state);
        let mut temp = state.clone();
        temp.axpy(dt.clone() * half.clone(), &k1, T::one());

        let k2 = f(t + dt.clone() * half, &temp);
        state.axpy(dt, &k2, T::one());

        Ok(())
    }

    fn order(&self) -> usize {
        2
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

/// Runge-Kutta 4th order (RK4) time integration
pub struct RungeKutta4;

impl<T: RealField + FromPrimitive> TimeIntegrator<T> for RungeKutta4 {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let two = T::from_f64(2.0).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation)
        })?;
        let six = T::from_f64(6.0).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation)
        })?;
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation)
        })?;
        
        let k1 = f(t.clone(), state);
        
        let mut temp = state.clone();
        temp.axpy(dt.clone() * half.clone(), &k1, T::one());
        let k2 = f(t.clone() + dt.clone() * half.clone(), &temp);
        
        temp.clone_from(state);
        temp.axpy(dt.clone() * half.clone(), &k2, T::one());
        let k3 = f(t.clone() + dt.clone() * half, &temp);
        
        temp.clone_from(state);
        temp.axpy(dt.clone(), &k3, T::one());
        let k4 = f(t + dt.clone(), &temp);
        
        // y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        state.axpy(dt.clone() / six.clone(), &k1, T::one());
        state.axpy(dt.clone() * two.clone() / six.clone(), &k2, T::one());
        state.axpy(dt.clone() * two / six.clone(), &k3, T::one());
        state.axpy(dt / six, &k4, T::one());
        
        Ok(())
    }

    fn order(&self) -> usize {
        4
    }

    fn is_explicit(&self) -> bool {
        true
    }
}

/// Backward Euler (implicit) time integration
pub struct BackwardEuler<T: RealField> {
    /// Tolerance for nonlinear solver
    pub tolerance: T,
    /// Maximum iterations for nonlinear solver
    pub max_iterations: usize,
}

impl<T: RealField + FromPrimitive> Default for BackwardEuler<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FromPrimitive> BackwardEuler<T> {
    /// Create a new BackwardEuler integrator with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the maximum number of iterations
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set the tolerance
    pub fn with_tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = tolerance;
        self
    }
}

impl<T: RealField> TimeIntegrator<T> for BackwardEuler<T> {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Implement implicit Backward Euler using fixed-point iteration
        // Solve: y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})

        let mut y_new = state.clone();
        let y_old = state.clone();
        let t_new = t + dt.clone();

        // Check for valid iteration count
        if self.max_iterations == 0 {
            return Err(crate::error::Error::InvalidConfiguration(
                "BackwardEuler requires max_iterations > 0".to_string()
            ));
        }

        // Fixed-point iteration: y_{n+1}^{k+1} = y_n + dt * f(t_{n+1}, y_{n+1}^k)
        for iteration in 0..self.max_iterations {
            let f_val = f(t_new.clone(), &y_new);
            let y_next = &y_old + &(&f_val * dt.clone());

            // Check convergence
            let error = (&y_next - &y_new).norm();
            if error < self.tolerance {
                *state = y_next;
                return Ok(());
            }

            y_new = y_next;

            // Prevent infinite loops
            if iteration == self.max_iterations - 1 {
                return Err(crate::error::Error::Convergence(
                    crate::error::ConvergenceErrorKind::MaxIterationsExceeded { 
                        max: self.max_iterations 
                    }
                ));
            }
        }

        Ok(())
    }

    fn order(&self) -> usize {
        1
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

/// Crank-Nicolson (implicit) time integration
pub struct CrankNicolson<T: RealField> {
    /// Tolerance for nonlinear solver
    pub tolerance: T,
    /// Maximum iterations for nonlinear solver
    pub max_iterations: usize,
}

impl<T: RealField + FromPrimitive> Default for CrankNicolson<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FromPrimitive> CrankNicolson<T> {
    /// Create a new CrankNicolson integrator with default settings
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the maximum number of iterations
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set the tolerance
    pub fn with_tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = tolerance;
        self
    }
}

impl<T: RealField> TimeIntegrator<T> for CrankNicolson<T> {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Implement Crank-Nicolson using fixed-point iteration
        // Solve: y_{n+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))

        let mut y_new = state.clone();
        let y_old = state.clone();
        let t_new = t.clone() + dt.clone();
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation)
        })?;
        let half_dt = dt * half;

        // Compute f(t_n, y_n) once
        let f_old = f(t, &y_old);

        // Check for valid iteration count
        if self.max_iterations == 0 {
            return Err(crate::error::Error::InvalidConfiguration(
                "CrankNicolson requires max_iterations > 0".to_string()
            ));
        }

        // Fixed-point iteration: y_{n+1}^{k+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}^k))
        for iteration in 0..self.max_iterations {
            let f_new = f(t_new.clone(), &y_new);
            let y_next = &y_old + &((&f_old + &f_new) * half_dt.clone());

            // Check convergence
            let error = (&y_next - &y_new).norm();
            if error < self.tolerance {
                *state = y_next;
                return Ok(());
            }

            y_new = y_next;

            // Prevent infinite loops
            if iteration == self.max_iterations - 1 {
                return Err(crate::error::Error::Convergence(
                    crate::error::ConvergenceErrorKind::MaxIterationsExceeded { 
                        max: self.max_iterations 
                    }
                ));
            }
        }

        Ok(())
    }

    fn order(&self) -> usize {
        2
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

/// Adaptive time stepping controller
/// Adaptive time step controller for error-based time step adjustment
pub struct AdaptiveTimeStepController<T: RealField> {
    /// Target error tolerance
    pub target_error: T,
    /// Safety factor for time step adjustment
    pub safety_factor: T,
    /// Maximum time step increase factor
    pub max_increase: T,
    /// Minimum time step decrease factor
    pub min_decrease: T,
}

impl<T: RealField + FromPrimitive> Default for AdaptiveTimeStepController<T> {
    fn default() -> Self {
        Self {
            target_error: T::from_f64(1e-6).unwrap_or_else(T::zero),
            safety_factor: T::from_f64(0.9).unwrap_or_else(T::zero),
            max_increase: T::from_f64(2.0).unwrap_or_else(T::zero),
            min_decrease: T::from_f64(0.1).unwrap_or_else(T::zero),
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> AdaptiveTimeStepController<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> T {
        let factor = (self.target_error / error).powf(T::one() / T::from_usize(order + 1).unwrap_or_else(T::one));
        let factor = factor * self.safety_factor;
        let factor = factor.min(self.max_increase).max(self.min_decrease);
        current_dt * factor
    }
}

pub struct VariableTimeStep<T: RealField> {
    /// Minimum allowed time step
    pub dt_min: T,
    /// Maximum allowed time step
    pub dt_max: T,
    /// Safety factor for step size adjustment
    pub safety_factor: T,
    /// Target error tolerance
    pub target_error: T,
}

impl<T: RealField + FromPrimitive> Default for VariableTimeStep<T> {
    fn default() -> Self {
        Self {
            dt_min: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            dt_max: T::from_f64(0.1).unwrap_or_else(|| T::one()),
            safety_factor: T::from_f64(0.9).unwrap_or_else(|| T::one()),
            target_error: T::from_f64(1e-6).unwrap_or_else(|| T::one()),
        }
    }
}

impl<T: RealField + FromPrimitive + Float> VariableTimeStep<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> T {
        if error < T::epsilon() {
            return num_traits::Float::min(self.dt_max.clone(), current_dt * T::from_f64(2.0).unwrap_or_else(|| T::one()));
        }
        
        let exponent = T::one() / T::from_f64(order as f64).unwrap_or_else(|| T::one());
        let factor = self.safety_factor.clone()
            * num_traits::Float::powf(self.target_error.clone() / error, exponent);
        
        let current_dt = current_dt * factor;
        let max_dt = num_traits::Float::max(current_dt, self.dt_min.clone());
        num_traits::Float::min(max_dt, self.dt_max.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Standard state for testing
    #[derive(Clone)]
    #[allow(dead_code)]
    struct TestState(f64);

    impl std::ops::AddAssign for TestState {
        fn add_assign(&mut self, rhs: Self) {
            self.0 += rhs.0;
        }
    }

    impl std::ops::Add for TestState {
        type Output = Self;
        fn add(self, rhs: Self) -> Self::Output {
            TestState(self.0 + rhs.0)
        }
    }

    impl std::ops::Mul<f64> for TestState {
        type Output = Self;
        fn mul(self, rhs: f64) -> Self::Output {
            TestState(self.0 * rhs)
        }
    }

    #[test]
    fn test_forward_euler() {
        use approx::assert_abs_diff_eq;

        let integrator = ForwardEuler;
        
        // Test with a simple ODE: dy/dt = -y
        // Solution: y(t) = y0 * exp(-t)
        let mut state = nalgebra::DVector::from_element(1, 1.0);
        let dt = 0.1;
        
        // Define derivative function
        let derivative = |_t: f64, y: &nalgebra::DVector<f64>| -> nalgebra::DVector<f64> {
            -y
        };
        
        // Take one step
        integrator.step(&mut state, 0.0, dt, derivative).expect("CRITICAL: Add proper error handling");
        
        // After one step: y ≈ y0 * (1 - dt) = 1.0 * (1 - 0.1) = 0.9
        assert_abs_diff_eq!(state[0], 0.9, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_time_step() {
        let controller = AdaptiveTimeStepController::<f64>::default();
        
        // Error is less than target
        let current_dt = controller.calculate_dt(0.01, 1e-8, 2);
        assert!(current_dt > 0.01);
        
        // Error is greater than target
        let current_dt = controller.calculate_dt(0.01, 1e-4, 2);
        assert!(current_dt < 0.01);
    }

    #[test]
    fn test_backward_euler_zero_iterations() {
        let integrator = BackwardEuler::new().with_max_iterations(0);
        let mut state = nalgebra::DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 0.1;

        let f = |_t: f64, y: &nalgebra::DVector<f64>| -y;

        // Should return an error for zero iterations
        let result = integrator.step(&mut state, t, dt, f);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), crate::error::Error::InvalidConfiguration(_)));
    }

    #[test]
    fn test_crank_nicolson_zero_iterations() {
        let integrator = CrankNicolson::new().with_max_iterations(0);
        let mut state = nalgebra::DVector::from_element(1, 1.0);
        let t = 0.0;
        let dt = 0.1;

        let f = |_t: f64, y: &nalgebra::DVector<f64>| -y;

        // Should return an error for zero iterations
        let result = integrator.step(&mut state, t, dt, f);
        assert!(result.is_err());
        assert!(matches!(result.unwrap_err(), crate::error::Error::InvalidConfiguration(_)));
    }
}