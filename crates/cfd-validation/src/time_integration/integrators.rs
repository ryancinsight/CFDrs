//! Time integration methods for validation.

use cfd_core::conversion::SafeFromF64;
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

impl<T: RealField + Copy + FromPrimitive + SafeFromF64> TimeIntegratorTrait<T> for RungeKutta2 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t, y);
        
        // Zero-copy: Reuse k1 vector for y_intermediate calculation
        let mut y_intermediate = k1.clone(); // Clone only k1, not y
        y_intermediate *= dt;
        y_intermediate += &*y;
        
        let k2 = f(t + dt, &y_intermediate);

        let half = T::from_f64_or_zero(HALF);
        *y += (k1 + k2) * (dt * half);
        Ok(())
    }

    fn order(&self) -> usize {
        2
    }
}

/// Fourth-order Runge-Kutta time integrator (RK4)
pub struct RungeKutta4;

impl<T: RealField + Copy + FromPrimitive + SafeFromF64> TimeIntegratorTrait<T> for RungeKutta4 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let half = T::from_f64_or_zero(HALF);
        
        let k1 = f(t, y);
        
        // Zero-copy: Reuse k1 for y_k1 calculation
        let mut y_k1 = k1.clone(); // Clone k1, not y
        y_k1 *= dt * half;
        y_k1 += &*y;
        
        let k2 = f(t + dt * half, &y_k1);
        
        // Zero-copy: Reuse k2 for y_k2 calculation
        let mut y_k2 = k2.clone(); // Clone k2, not y
        y_k2 *= dt * half;
        y_k2 += &*y;
        
        let k3 = f(t + dt * half, &y_k2);
        
        // Zero-copy: Reuse k3 for y_k3 calculation
        let mut y_k3 = k3.clone(); // Clone k3, not y
        y_k3 *= dt;
        y_k3 += &*y;
        
        let k4 = f(t + dt, &y_k3);

        let sixth = T::from_f64_or_zero(ONE_SIXTH);
        let two = T::from_f64_or_zero(TWO);
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
