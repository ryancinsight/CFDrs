//! Time integration schemes.

use crate::Result;
use nalgebra::RealField;

/// Trait for time integration schemes
pub trait TimeIntegrator<T: RealField>: Send + Sync {
    /// State type
    type State;

    /// Perform one time step
    fn step(
        &self,
        state: &mut Self::State,
        dt: T,
        derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()>;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Check if the scheme is explicit
    fn is_explicit(&self) -> bool;
}

/// Forward Euler (explicit) time integration
pub struct ForwardEuler;

impl<T: RealField> TimeIntegrator<T> for ForwardEuler {
    type State = nalgebra::DVector<T>;

    fn step(
        &self,
        state: &mut Self::State,
        dt: T,
        derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()> {
        let k1 = derivative(state)?;
        state.axpy(dt, &k1, T::one());
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

impl<T: RealField> TimeIntegrator<T> for RungeKutta2 
where
    T: From<f64>,
{
    type State = nalgebra::DVector<T>;

    fn step(
        &self,
        state: &mut Self::State,
        dt: T,
        derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()> {
        let k1 = derivative(state)?;
        let mut temp = state.clone();
        temp.axpy(dt * T::from(0.5).unwrap(), &k1, T::one());
        let k2 = derivative(&temp)?;
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

impl<T: RealField> TimeIntegrator<T> for RungeKutta4
where
    T: From<f64>,
{
    type State = nalgebra::DVector<T>;

    fn step(
        &self,
        state: &mut Self::State,
        dt: T,
        derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()> {
        let two = T::from(2.0).unwrap();
        let six = T::from(6.0).unwrap();
        let half = T::from(0.5).unwrap();
        
        let k1 = derivative(state)?;
        let mut temp1 = state.clone();
        temp1.axpy(dt * half, &k1, T::one());
        
        let k2 = derivative(&temp1)?;
        let mut temp2 = state.clone();
        temp2.axpy(dt * half, &k2, T::one());
        
        let k3 = derivative(&temp2)?;
        let mut temp3 = state.clone();
        temp3.axpy(dt, &k3, T::one());
        
        let k4 = derivative(&temp3)?;
        
        // state += (k1 + 2*k2 + 2*k3 + k4) * dt/6
        let dt_over_six = dt / six;
        state.axpy(dt_over_six, &k1, T::one());
        state.axpy(dt_over_six * two, &k2, T::one());
        state.axpy(dt_over_six * two, &k3, T::one());
        state.axpy(dt_over_six, &k4, T::one());
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

impl<T: RealField> Default for BackwardEuler<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from(1e-10).unwrap(),
            max_iterations: 100,
        }
    }
}

impl<T: RealField> TimeIntegrator<T> for BackwardEuler<T> {
    type State = nalgebra::DVector<T>;

    fn step(
        &self,
        _state: &mut Self::State,
        _dt: T,
        _derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()> {
        // Note: This is a placeholder implementation
        // Actual implementation would require a nonlinear solver
        // to solve: state_new = state_old + dt * f(state_new)
        tracing::warn!("Backward Euler requires problem-specific implementation");
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

impl<T: RealField> Default for CrankNicolson<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from(1e-10).unwrap(),
            max_iterations: 100,
        }
    }
}

impl<T: RealField> TimeIntegrator<T> for CrankNicolson<T>
where
    T: From<f64>,
{
    type State = nalgebra::DVector<T>;

    fn step(
        &self,
        _state: &mut Self::State,
        _dt: T,
        _derivative: impl Fn(&Self::State) -> Result<Self::State>,
    ) -> Result<()> {
        // Note: This is a placeholder implementation
        // Actual implementation would require a nonlinear solver
        // to solve: state_new = state_old + dt/2 * (f(state_old) + f(state_new))
        tracing::warn!("Crank-Nicolson requires problem-specific implementation");
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
pub struct AdaptiveTimeStep<T: RealField> {
    /// Minimum time step
    pub dt_min: T,
    /// Maximum time step
    pub dt_max: T,
    /// Safety factor
    pub safety_factor: T,
    /// Target relative error
    pub target_error: T,
}

impl<T: RealField> Default for AdaptiveTimeStep<T> {
    fn default() -> Self {
        Self {
            dt_min: T::from(1e-10).unwrap(),
            dt_max: T::from(0.1).unwrap(),
            safety_factor: T::from(0.9).unwrap(),
            target_error: T::from(1e-6).unwrap(),
        }
    }
}

impl<T: RealField> AdaptiveTimeStep<T> {
    /// Compute new time step based on error estimate
    pub fn compute_dt(&self, current_dt: T, error: T, order: usize) -> T {
        if error < T::epsilon() {
            return (self.dt_max).min(current_dt * T::from(2.0).unwrap());
        }
        
        let exponent = T::one() / T::from(order as f64).unwrap();
        let factor = self.safety_factor 
            * num_traits::Float::powf(self.target_error / error, exponent);
        
        let new_dt = current_dt * factor;
        new_dt.max(self.dt_min).min(self.dt_max)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Simple state for testing
    #[derive(Clone)]
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
        let integrator = ForwardEuler;
        let mut state = TestState(1.0);
        
        // dy/dt = -y
        let derivative = |s: &TestState| Ok(TestState(-s.0));
        
        integrator.step(&mut state, 0.1, derivative).unwrap();
        assert!((state.0 - 0.9).abs() < 1e-10);
    }

    #[test]
    fn test_adaptive_time_step() {
        let controller = AdaptiveTimeStep::<f64>::default();
        
        // Error is less than target
        let new_dt = controller.compute_dt(0.01, 1e-8, 2);
        assert!(new_dt > 0.01);
        
        // Error is greater than target
        let new_dt = controller.compute_dt(0.01, 1e-4, 2);
        assert!(new_dt < 0.01);
    }
}