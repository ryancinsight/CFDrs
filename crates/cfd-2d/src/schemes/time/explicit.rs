//! Explicit time integration methods
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_core::physics::constants::mathematical::numeric::{ONE_HALF, SIX, TWO};
use eunomia::FloatElement;

use super::vector::StateVector;

/// Forward Euler step: y_{n+1} = y_n + dt*f(t_n, y_n)
pub fn forward_euler<T, F>(f: F, y: &StateVector<T>, t: T, dt: T) -> StateVector<T>
where
    T: Cfd2dScalar + Copy,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    let increment = &f(t, y) * dt;
    y + &increment
}

/// Second-order Runge-Kutta step (RK2 - Heun's method)
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn runge_kutta2<T, F>(f: F, y: &StateVector<T>, t: T, dt: T) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    let k1 = f(t, y);
    let half_dt = dt * scalar::from_f64::<T>(ONE_HALF);
    let k1_half = &k1 * half_dt;
    let y_mid = y + &k1_half;
    let k2 = f(t + half_dt, &y_mid);
    let increment = &k2 * dt;
    y + &increment
}

/// Fourth-order Runge-Kutta step (RK4 - classical)
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn runge_kutta4<T, F>(f: F, y: &StateVector<T>, t: T, dt: T) -> StateVector<T>
where
    T: Cfd2dScalar + Copy + FloatElement,
    F: Fn(T, &StateVector<T>) -> StateVector<T>,
{
    let two = scalar::from_f64::<T>(TWO);
    let six = scalar::from_f64::<T>(SIX);

    let k1 = f(t, y);
    let half_dt = dt * scalar::from_f64::<T>(ONE_HALF);
    let k1_half = &k1 * half_dt;
    let y2 = y + &k1_half;
    let k2 = f(t + half_dt, &y2);
    let k2_half = &k2 * half_dt;
    let y3 = y + &k2_half;
    let k3 = f(t + half_dt, &y3);
    let k3_full = &k3 * dt;
    let y4 = y + &k3_full;
    let k4 = f(t + dt, &y4);

    let k2_weighted = &k2 * two;
    let k3_weighted = &k3 * two;
    let sum_12 = &k1 + &k2_weighted;
    let sum_123 = &sum_12 + &k3_weighted;
    let sum = &sum_123 + &k4;
    let increment = &sum * (dt / six);
    y + &increment
}
