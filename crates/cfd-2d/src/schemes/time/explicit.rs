//! Explicit time integration methods

use cfd_core::constants::mathematical::numeric::{ONE_HALF, SIX, TWO};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// Forward Euler step: y_{n+1} = y_n + dt*f(t_n, y_n)
pub fn forward_euler<T, F>(f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    y + f(t, y) * dt
}

/// Second-order Runge-Kutta step (RK2 - Heun's method)
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn runge_kutta2<T, F>(f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    let k1 = f(t, y);
    let half_dt = dt * T::from_f64(ONE_HALF).unwrap_or_else(T::zero);
    let y_mid = y + &k1 * half_dt;
    let k2 = f(t + half_dt, &y_mid);
    y + k2 * dt
}

/// Fourth-order Runge-Kutta step (RK4 - classical)
///
/// Reference: Butcher (2016) - Numerical Methods for Ordinary Differential Equations
pub fn runge_kutta4<T, F>(f: F, y: &DVector<T>, t: T, dt: T) -> DVector<T>
where
    T: RealField + Copy + FromPrimitive,
    F: Fn(T, &DVector<T>) -> DVector<T>,
{
    let two = T::from_f64(TWO).unwrap_or_else(T::zero);
    let six = T::from_f64(SIX).unwrap_or_else(T::zero);

    let k1 = f(t, y);
    let half_dt = dt * T::from_f64(ONE_HALF).unwrap_or_else(T::zero);
    let y2 = y + &k1 * half_dt;
    let k2 = f(t + half_dt, &y2);
    let y3 = y + &k2 * half_dt;
    let k3 = f(t + half_dt, &y3);
    let y4 = y + &k3 * dt;
    let k4 = f(t + dt, &y4);

    y + (k1 + k2 * two + k3 * two + k4) * (dt / six)
}
