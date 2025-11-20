//! Time integration schemes
//!
//! This module provides various time integration methods for solving ordinary differential equations (ODEs).
//! The schemes are organized by type:
//! - Explicit methods: Forward Euler, Runge-Kutta (RK2, RK4)
//! - Implicit methods: Backward Euler, Crank-Nicolson
//! - Multi-step methods: Adams-Bashforth 2, BDF2, BDF3

use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

mod adaptive;
mod explicit;
mod implicit;
mod multistep;
#[cfg(test)]
mod tests;
mod types;

pub use adaptive::{AdaptationStrategy, AdaptiveController, AdaptiveTimeIntegrator};
pub use types::TimeScheme;

/// Time integrator for ODEs
pub struct TimeIntegrator<T: RealField + Copy> {
    scheme: TimeScheme,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Clone> TimeIntegrator<T> {
    /// Create new time integrator
    #[must_use]
    pub fn new(scheme: TimeScheme) -> Self {
        Self {
            scheme,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Perform time step
    ///
    /// # Implementation Notes
    ///
    /// - **Implicit schemes** (BackwardEuler, CrankNicolson): Use fixed-point iteration
    ///   to solve the implicit equations. For stiff systems, consider using
    ///   `step_with_history()` with BDF2/BDF3 for better stability.
    ///
    /// - **Multi-step schemes** (AdamsBashforth2): Require history. Use
    ///   `step_with_history()` with previous solution for proper implementation.
    ///   This method falls back to Forward Euler for first step.
    ///
    /// # Arguments
    ///
    /// * `f` - Right-hand side function f(t, y) for dy/dt = f(t, y)
    /// * `y` - Current solution at time t
    /// * `t` - Current time
    /// * `dt` - Time step size
    ///
    /// # Returns
    ///
    /// Solution at time t + dt
    pub fn step<F>(&self, f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self.scheme {
            TimeScheme::ForwardEuler => explicit::forward_euler(f, y, t, dt),
            TimeScheme::BackwardEuler => implicit::backward_euler(f, y, t, dt),
            TimeScheme::CrankNicolson => implicit::crank_nicolson(f, y, t, dt),
            TimeScheme::RungeKutta2 => explicit::runge_kutta2(f, y, t, dt),
            TimeScheme::RungeKutta4 => explicit::runge_kutta4(f, y, t, dt),
            TimeScheme::AdamsBashforth2 => {
                // Adams-Bashforth 2nd order requires history.
                // For single-step call, fall back to RK2.
                // Use step_with_history() for proper multi-step implementation.
                // Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
                explicit::runge_kutta2(f, y, t, dt)
            }
            TimeScheme::BDF2 | TimeScheme::BDF3 => {
                // BDF schemes require history - delegate to step_with_history()
                // For single-step call, fall back to Backward Euler
                self.step_with_history(f, y, None, None, t, dt)
            }
        }
    }

    /// Get scheme order of accuracy
    #[must_use]
    pub fn order(&self) -> usize {
        self.scheme.order()
    }

    /// Check if scheme is explicit
    #[must_use]
    pub fn is_explicit(&self) -> bool {
        self.scheme.is_explicit()
    }

    /// Perform time step with history (required for multi-step methods)
    ///
    /// # Arguments
    ///
    /// * `f` - Right-hand side function f(t, y) for dy/dt = f(t, y)
    /// * `y_curr` - Current solution at time t_n
    /// * `y_prev` - Previous solution at time t_{n-1} (required for multi-step methods)
    /// * `y_prev2` - Second previous solution at time t_{n-2} (required for BDF3)
    /// * `t` - Current time t_n
    /// * `dt` - Time step size
    ///
    /// # Returns
    ///
    /// Solution at time t_{n+1}
    ///
    /// # Implementation Notes
    ///
    /// - **BDF2/BDF3**: A-stable (BDF2) or stiffly-stable (BDF3) implicit methods
    ///   suitable for stiff systems. If history not provided, falls back to
    ///   Backward Euler.
    ///
    /// - **Adams-Bashforth2**: 2nd-order explicit multi-step method. If history
    ///   not provided, falls back to RK2.
    ///
    /// - **Other schemes**: History is ignored, delegates to `step()`.
    ///
    /// # References
    ///
    /// - Curtiss & Hirschfelder (1952): BDF formulas
    /// - Butcher (2016): Adams-Bashforth methods
    /// - Hairer & Wanner (1996): Stiff ODEs
    pub fn step_with_history<F>(
        &self,
        f: F,
        y_curr: &DVector<T>,
        y_prev: Option<&DVector<T>>,
        y_prev2: Option<&DVector<T>>,
        t: T,
        dt: T,
    ) -> DVector<T>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self.scheme {
            TimeScheme::AdamsBashforth2 => multistep::adams_bashforth2(f, y_curr, y_prev, t, dt),
            TimeScheme::BDF2 => multistep::bdf2(f, y_curr, y_prev, t, dt),
            TimeScheme::BDF3 => multistep::bdf3(f, y_curr, y_prev, y_prev2, t, dt),
            _ => {
                // For other schemes, delegate to regular step() method
                self.step(f, y_curr, t, dt)
            }
        }
    }
}
