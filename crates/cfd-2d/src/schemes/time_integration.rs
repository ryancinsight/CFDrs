//! Time integration schemes

use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Time integration scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TimeScheme {
    /// Forward Euler (explicit)
    ForwardEuler,
    /// Backward Euler (implicit)
    BackwardEuler,
    /// Crank-Nicolson
    CrankNicolson,
    /// Second-order Runge-Kutta
    RungeKutta2,
    /// Fourth-order Runge-Kutta
    RungeKutta4,
    /// Adams-Bashforth (2nd order)
    AdamsBashforth2,
}

/// Time integrator for ODEs
pub struct TimeIntegrator<T: RealField + Copy> {
    scheme: TimeScheme,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive + Clone> TimeIntegrator<T> {
    /// Create new time integrator
    pub fn new(scheme: TimeScheme) -> Self {
        Self {
            scheme,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Perform time step
    pub fn step<F>(&self, f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self.scheme {
            TimeScheme::ForwardEuler => {
                y + f(t, y) * dt
            },
            TimeScheme::RungeKutta2 => {
                let k1 = f(t, y);
                let half_dt = dt / T::from_f64(2.0).unwrap_or_else(T::zero);
                let y_mid = y + &k1 * half_dt;
                let k2 = f(t + half_dt, &y_mid);
                y + k2 * dt
            },
            TimeScheme::RungeKutta4 => {
                let two = T::from_f64(2.0).unwrap_or_else(T::zero);
                let six = T::from_f64(6.0).unwrap_or_else(T::zero);
                
                let k1 = f(t, y);
                let half_dt = dt / two;
                let y2 = y + &k1 * half_dt;
                let k2 = f(t + half_dt, &y2);
                let y3 = y + &k2 * half_dt;
                let k3 = f(t + half_dt, &y3);
                let y4 = y + &k3 * dt;
                let k4 = f(t + dt, &y4);
                
                y + (k1 + k2 * two + k3 * two + k4) * (dt / six)
            },
            _ => {
                // Default to Forward Euler for unimplemented schemes
                y + f(t, y) * dt
            }
        }
    }
    
    /// Get scheme order of accuracy
    pub fn order(&self) -> usize {
        match self.scheme {
            TimeScheme::ForwardEuler | TimeScheme::BackwardEuler => 1,
            TimeScheme::CrankNicolson | TimeScheme::RungeKutta2 | TimeScheme::AdamsBashforth2 => 2,
            TimeScheme::RungeKutta4 => 4,
        }
    }
    
    /// Check if scheme is explicit
    pub fn is_explicit(&self) -> bool {
        matches!(
            self.scheme,
            TimeScheme::ForwardEuler | TimeScheme::RungeKutta2 | TimeScheme::RungeKutta4 | TimeScheme::AdamsBashforth2
        )
    }
}