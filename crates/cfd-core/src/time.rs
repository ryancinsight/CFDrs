//! Time integration schemes.

use crate::Result;
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
        let k1 = f(t.clone(), state);
        let mut temp = state.clone();
        temp.axpy(dt.clone() * T::from_f64(0.5).unwrap(), &k1, T::one());
        
        let k2 = f(t + dt.clone() * T::from_f64(0.5).unwrap(), &temp);
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
        let two = T::from_f64(2.0).unwrap();
        let six = T::from_f64(6.0).unwrap();
        let half = T::from_f64(0.5).unwrap();
        
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
            tolerance: T::from_f64(1e-10).unwrap(),
            max_iterations: 100,
        }
    }
}

impl<T: RealField> TimeIntegrator<T> for BackwardEuler<T> {
    type State = DVector<T>;

    fn step<F>(&self, _state: &mut Self::State, _t: T, _dt: T, _f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // TODO: Implement implicit solver
        // This requires solving: y_{n+1} = y_n + dt * f(t_{n+1}, y_{n+1})
        Err(crate::Error::NotImplemented(
            "Backward Euler requires nonlinear solver".to_string(),
        ))
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
            tolerance: T::from_f64(1e-10).unwrap(),
            max_iterations: 100,
        }
    }
}

impl<T: RealField> TimeIntegrator<T> for CrankNicolson<T> {
    type State = DVector<T>;

    fn step<F>(&self, _state: &mut Self::State, _t: T, _dt: T, _f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // TODO: Implement implicit solver
        // This requires solving: y_{n+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))
        Err(crate::Error::NotImplemented(
            "Crank-Nicolson requires nonlinear solver".to_string(),
        ))
    }

    fn order(&self) -> usize {
        2
    }

    fn is_explicit(&self) -> bool {
        false
    }
}

/// Adaptive time stepping controller
pub struct AdaptiveTimeStep<T: RealField> {
    /// Minimum allowed time step
    pub dt_min: T,
    /// Maximum allowed time step
    pub dt_max: T,
    /// Safety factor for step size adjustment
    pub safety_factor: T,
    /// Target error tolerance
    pub target_error: T,
}

impl<T: RealField + FromPrimitive> Default for AdaptiveTimeStep<T> {
    fn default() -> Self {
        Self {
            dt_min: T::from_f64(1e-10).unwrap(),
            dt_max: T::from_f64(0.1).unwrap(),
            safety_factor: T::from_f64(0.9).unwrap(),
            target_error: T::from_f64(1e-6).unwrap(),
        }
    }
}

impl<T: RealField + FromPrimitive + Float> AdaptiveTimeStep<T> {
    /// Calculate new time step based on error estimate
    pub fn calculate_dt(&self, current_dt: T, error: T, order: usize) -> T {
        if error < T::epsilon() {
            return num_traits::Float::min(self.dt_max.clone(), current_dt * T::from_f64(2.0).unwrap());
        }
        
        let exponent = T::one() / T::from_f64(order as f64).unwrap();
        let factor = self.safety_factor.clone()
            * num_traits::Float::powf(self.target_error.clone() / error, exponent);
        
        let new_dt = current_dt * factor;
        let max_dt = num_traits::Float::max(new_dt, self.dt_min.clone());
        num_traits::Float::min(max_dt, self.dt_max.clone())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Simple state for testing
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
        integrator.step(&mut state, 0.0, dt, derivative).unwrap();
        
        // After one step: y ≈ y0 * (1 - dt) = 1.0 * (1 - 0.1) = 0.9
        assert_abs_diff_eq!(state[0], 0.9, epsilon = 1e-10);
    }

    #[test]
    fn test_adaptive_time_step() {
        let controller = AdaptiveTimeStep::<f64>::default();
        
        // Error is less than target
        let new_dt = controller.calculate_dt(0.01, 1e-8, 2);
        assert!(new_dt > 0.01);
        
        // Error is greater than target
        let new_dt = controller.calculate_dt(0.01, 1e-4, 2);
        assert!(new_dt < 0.01);
    }
}