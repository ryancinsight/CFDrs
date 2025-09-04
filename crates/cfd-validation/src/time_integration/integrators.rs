//! Time integration methods for validation.

use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;
use std::marker::PhantomData;

// Integration constants from numerical analysis literature
/// One half constant for improved numerical stability
pub const HALF: f64 = 0.5;
/// One third constant for Runge-Kutta calculations
pub const ONE_THIRD: f64 = 1.0 / 3.0;
/// One sixth constant for Simpson's rule and RK4
pub const ONE_SIXTH: f64 = 1.0 / 6.0;
/// Two constant for numerical calculations
pub const TWO: f64 = 2.0;
/// Four constant for fourth-order methods
pub const FOUR: f64 = 4.0;

/// Time integrator trait for validation
pub trait TimeIntegratorTrait<T: RealField + Copy> {
    /// Take one time step
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>;

    /// Get the order of accuracy
    fn order(&self) -> usize;
}

/// Forward Euler time integrator (first-order explicit method)
pub struct ForwardEuler;

impl<T: RealField + Copy + FromPrimitive> TimeIntegratorTrait<T> for ForwardEuler {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t, y);
        *y += k1 * dt;
        Ok(())
    }

    fn order(&self) -> usize {
        1
    }
}

/// Second-order Runge-Kutta time integrator (Heun's method)
pub struct RungeKutta2;

impl<T: RealField + Copy + FromPrimitive> TimeIntegratorTrait<T> for RungeKutta2 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t, y);
        let y_intermediate = y.clone() + k1.clone() * dt;
        let k2 = f(t + dt, &y_intermediate);

        let half = T::from_f64(HALF).unwrap_or_else(|| T::zero());
        *y += (k1 + k2) * (dt * half);
        Ok(())
    }

    fn order(&self) -> usize {
        2
    }
}

/// Fourth-order Runge-Kutta time integrator (RK4)
pub struct RungeKutta4;

impl<T: RealField + Copy + FromPrimitive> TimeIntegratorTrait<T> for RungeKutta4 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let half = T::from_f64(HALF).unwrap_or_else(|| T::zero());
        let k1 = f(t, y);
        let y_k1 = y.clone() + k1.clone() * (dt * half);
        let k2 = f(t + dt * half, &y_k1);
        let y_k2 = y.clone() + k2.clone() * (dt * half);
        let k3 = f(t + dt * half, &y_k2);
        let y_k3 = y.clone() + k3.clone() * dt;
        let k4 = f(t + dt, &y_k3);

        let sixth = T::from_f64(ONE_SIXTH).unwrap_or_else(|| T::zero());
        let two = T::from_f64(TWO).unwrap_or_else(|| T::zero());
        *y += (k1 + k2 * two + k3 * two + k4) * (dt * sixth);
        Ok(())
    }

    fn order(&self) -> usize {
        4
    }
}

/// Enum wrapper for time integrators
pub enum TimeIntegratorEnum<T: RealField + Copy> {
    /// Forward Euler
    ForwardEuler(ForwardEuler, PhantomData<T>),
    /// Second-order Runge-Kutta
    RungeKutta2(RungeKutta2, PhantomData<T>),
    /// Fourth-order Runge-Kutta
    RungeKutta4(RungeKutta4, PhantomData<T>),
}

impl<T: RealField + Copy + FromPrimitive> TimeIntegratorTrait<T> for TimeIntegratorEnum<T> {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self {
            TimeIntegratorEnum::ForwardEuler(integrator, _) => integrator.step(y, t, dt, f),
            TimeIntegratorEnum::RungeKutta2(integrator, _) => integrator.step(y, t, dt, f),
            TimeIntegratorEnum::RungeKutta4(integrator, _) => integrator.step(y, t, dt, f),
        }
    }

    fn order(&self) -> usize {
        match self {
            TimeIntegratorEnum::ForwardEuler(_, _) => 1,
            TimeIntegratorEnum::RungeKutta2(_, _) => 2,
            TimeIntegratorEnum::RungeKutta4(_, _) => 4,
        }
    }
}
