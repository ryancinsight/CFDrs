//! Time integration schemes.

use crate::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::cast::FromPrimitive;

/// Trait for time integration schemes
pub trait TimeIntegrator<T: RealField + Copy>: Send + Sync {
    /// State type
    type State;

    /// Perform one time step
    ///
    /// # Errors
    /// Returns an error if the scheme fails its internal numerical checks.
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

impl<T: RealField + Copy> TimeIntegrator<T> for ForwardEuler {
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

impl<T: RealField + FromPrimitive + Copy> TimeIntegrator<T> for RungeKutta2 {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidValue {
                value: "Failed to convert 0.5 to T".to_string(),
            })
        })?;

        let k1 = f(t, state);
        let mut state_buffer = state.clone();
        state_buffer.axpy(dt * half, &k1, T::one());

        let k2 = f(t + dt * half, &state_buffer);
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

impl<T: RealField + FromPrimitive + Copy> TimeIntegrator<T> for RungeKutta4 {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let two = T::from_f64(2.0).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidValue {
                value: "Failed to convert 2.0 to T".to_string(),
            })
        })?;
        let six = T::from_f64(6.0).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidValue {
                value: "Failed to convert 6.0 to T".to_string(),
            })
        })?;
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidValue {
                value: "Failed to convert 0.5 to T".to_string(),
            })
        })?;

        let k1 = f(t, state);

        let mut state_buffer = state.clone();
        state_buffer.axpy(dt * half, &k1, T::one());
        let k2 = f(t + dt * half, &state_buffer);

        state_buffer.clone_from(state);
        state_buffer.axpy(dt * half, &k2, T::one());
        let k3 = f(t + dt * half, &state_buffer);

        state_buffer.clone_from(state);
        state_buffer.axpy(dt, &k3, T::one());
        let k4 = f(t + dt, &state_buffer);

        // y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        state.axpy(dt / six, &k1, T::one());
        state.axpy(dt * two / six, &k2, T::one());
        state.axpy(dt * two / six, &k3, T::one());
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
pub struct BackwardEuler<T: RealField + Copy> {
    /// Tolerance for nonlinear solver
    pub tolerance: T,
    /// Maximum iterations for nonlinear solver
    pub max_iterations: usize,
}

impl<T: RealField + FromPrimitive + Copy> Default for BackwardEuler<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> BackwardEuler<T> {
    /// Create a `BackwardEuler` integrator with default settings
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the maximum number of iterations
    #[must_use]
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set the tolerance
    #[must_use]
    pub fn with_tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = tolerance;
        self
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

        let mut next_candidate = state.clone();
        let previous_state = state.clone();
        let t_next = t + dt;

        // Check for valid iteration count
        if self.max_iterations == 0 {
            return Err(crate::error::Error::InvalidConfiguration(
                "BackwardEuler requires max_iterations > 0".to_string(),
            ));
        }

        // Fixed-point iteration: y_{n+1}^{k+1} = y_n + dt * f(t_{n+1}, y_{n+1}^k)
        for iteration in 0..self.max_iterations {
            let f_val = f(t_next, &next_candidate);
            let next_state = previous_state.clone() + f_val * dt;

            // Check convergence
            let error = (&next_state - &next_candidate).norm();
            if error < self.tolerance {
                *state = next_state;
                return Ok(());
            }

            next_candidate = next_state;

            // Prevent infinite loops
            if iteration == self.max_iterations - 1 {
                return Err(crate::error::Error::Convergence(
                    crate::error::ConvergenceErrorKind::MaxIterationsExceeded {
                        max: self.max_iterations,
                    },
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
pub struct CrankNicolson<T: RealField + Copy> {
    /// Tolerance for nonlinear solver
    pub tolerance: T,
    /// Maximum iterations for nonlinear solver
    pub max_iterations: usize,
}

impl<T: RealField + FromPrimitive + Copy> Default for CrankNicolson<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from_f64(1e-10).unwrap_or_else(|| T::one()),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> CrankNicolson<T> {
    /// Create a `CrankNicolson` integrator with default settings
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the maximum number of iterations
    #[must_use]
    pub fn with_max_iterations(mut self, max_iterations: usize) -> Self {
        self.max_iterations = max_iterations;
        self
    }

    /// Set the tolerance
    #[must_use]
    pub fn with_tolerance(mut self, tolerance: T) -> Self {
        self.tolerance = tolerance;
        self
    }
}

impl<T: RealField + Copy> TimeIntegrator<T> for CrankNicolson<T> {
    type State = DVector<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Implement Crank-Nicolson using fixed-point iteration
        // Solve: y_{n+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))

        let mut next_candidate = state.clone();
        let previous_state = state.clone();
        let t_next = t + dt;
        let half = T::from_f64(0.5).ok_or_else(|| {
            crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidValue {
                value: "Failed to convert 0.5 to T".to_string(),
            })
        })?;
        let half_dt = dt * half;

        // Compute f(t_n, y_n) once
        let f_previous = f(t, &previous_state);

        // Check for valid iteration count
        if self.max_iterations == 0 {
            return Err(crate::error::Error::InvalidConfiguration(
                "CrankNicolson requires max_iterations > 0".to_string(),
            ));
        }

        // Fixed-point iteration: y_{n+1}^{k+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}^k))
        for iteration in 0..self.max_iterations {
            let f_current = f(t_next, &next_candidate);
            let next_state = previous_state.clone() + (f_previous.clone() + f_current) * half_dt;

            // Check convergence
            let error = (&next_state - &next_candidate).norm();
            if error < self.tolerance {
                *state = next_state;
                return Ok(());
            }

            next_candidate = next_state;

            // Prevent infinite loops
            if iteration == self.max_iterations - 1 {
                return Err(crate::error::Error::Convergence(
                    crate::error::ConvergenceErrorKind::MaxIterationsExceeded {
                        max: self.max_iterations,
                    },
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
