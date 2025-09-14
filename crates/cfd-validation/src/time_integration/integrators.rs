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
    /// Take one time step with workspace buffer for zero-copy operation
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>;

    /// Take one time step with pre-allocated workspace for zero-copy performance
    fn step_with_workspace<F>(&self, y: &mut DVector<T>, workspace: &mut [DVector<T>], t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Get required workspace size for zero-copy operation
    fn workspace_size(&self) -> usize;
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

    fn step_with_workspace<F>(&self, y: &mut DVector<T>, _workspace: &mut [DVector<T>], t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        // Forward Euler doesn't need workspace - just delegate to regular step
        self.step(y, t, dt, f)
    }

    fn order(&self) -> usize {
        1
    }

    fn workspace_size(&self) -> usize {
        0 // Forward Euler requires no workspace
    }
}

/// Second-order Runge-Kutta time integrator (Heun's method)
pub struct RungeKutta2;

impl<T: RealField + Copy + FromPrimitive + SafeFromF64> TimeIntegratorTrait<T> for RungeKutta2 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        // Default implementation using internal workspace allocation
        let mut workspace = vec![DVector::zeros(y.len()); <Self as TimeIntegratorTrait<T>>::workspace_size(self)];
        self.step_with_workspace(y, &mut workspace, t, dt, f)
    }

    fn step_with_workspace<F>(&self, y: &mut DVector<T>, workspace: &mut [DVector<T>], t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        if workspace.len() < <Self as TimeIntegratorTrait<T>>::workspace_size(self) {
            return Err(cfd_core::error::Error::InvalidInput("Insufficient workspace size".to_string()));
        }

        let k1 = f(t, y);
        
        // True zero-copy: Use pre-allocated workspace buffer
        let y_intermediate = &mut workspace[0];
        y_intermediate.copy_from(&k1);
        *y_intermediate *= dt;
        *y_intermediate += &*y;
        
        let k2 = f(t + dt, y_intermediate);

        let half = T::from_f64_or_zero(HALF);
        *y += (k1 + k2) * (dt * half);
        Ok(())
    }

    fn order(&self) -> usize {
        2
    }

    fn workspace_size(&self) -> usize {
        1 // One intermediate vector for y_intermediate
    }
}

/// Fourth-order Runge-Kutta time integrator (RK4)
pub struct RungeKutta4;

impl<T: RealField + Copy + FromPrimitive + SafeFromF64> TimeIntegratorTrait<T> for RungeKutta4 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        // Default implementation using internal workspace allocation
        let mut workspace = vec![DVector::zeros(y.len()); <Self as TimeIntegratorTrait<T>>::workspace_size(self)];
        self.step_with_workspace(y, &mut workspace, t, dt, f)
    }

    fn step_with_workspace<F>(&self, y: &mut DVector<T>, workspace: &mut [DVector<T>], t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        if workspace.len() < <Self as TimeIntegratorTrait<T>>::workspace_size(self) {
            return Err(cfd_core::error::Error::InvalidInput("Insufficient workspace size".to_string()));
        }

        let half = T::from_f64_or_zero(HALF);
        
        let k1 = f(t, y);
        
        // True zero-copy: Use pre-allocated workspace buffers
        let y_k1 = &mut workspace[0];
        y_k1.copy_from(&k1);
        *y_k1 *= dt * half;
        *y_k1 += &*y;
        
        let k2 = f(t + dt * half, y_k1);
        
        let y_k2 = &mut workspace[1];
        y_k2.copy_from(&k2);
        *y_k2 *= dt * half;
        *y_k2 += &*y;
        
        let k3 = f(t + dt * half, y_k2);
        
        let y_k3 = &mut workspace[2];
        y_k3.copy_from(&k3);
        *y_k3 *= dt;
        *y_k3 += &*y;
        
        let k4 = f(t + dt, y_k3);

        let sixth = T::from_f64_or_zero(ONE_SIXTH);
        let two = T::from_f64_or_zero(TWO);
        *y += (k1 + k2 * two + k3 * two + k4) * (dt * sixth);
        Ok(())
    }

    fn order(&self) -> usize {
        4
    }

    fn workspace_size(&self) -> usize {
        3 // Three intermediate vectors: y_k1, y_k2, y_k3
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

impl<T: RealField + Copy + FromPrimitive + SafeFromF64> TimeIntegratorTrait<T> for TimeIntegratorEnum<T> {
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

    fn step_with_workspace<F>(&self, y: &mut DVector<T>, workspace: &mut [DVector<T>], t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        match self {
            TimeIntegratorEnum::ForwardEuler(integrator, _) => integrator.step_with_workspace(y, workspace, t, dt, f),
            TimeIntegratorEnum::RungeKutta2(integrator, _) => integrator.step_with_workspace(y, workspace, t, dt, f),
            TimeIntegratorEnum::RungeKutta4(integrator, _) => integrator.step_with_workspace(y, workspace, t, dt, f),
        }
    }

    fn order(&self) -> usize {
        match self {
            TimeIntegratorEnum::ForwardEuler(_, _) => 1,
            TimeIntegratorEnum::RungeKutta2(_, _) => 2,
            TimeIntegratorEnum::RungeKutta4(_, _) => 4,
        }
    }

    fn workspace_size(&self) -> usize {
        match self {
            TimeIntegratorEnum::ForwardEuler(integrator, _) => <ForwardEuler as TimeIntegratorTrait<T>>::workspace_size(integrator),
            TimeIntegratorEnum::RungeKutta2(integrator, _) => <RungeKutta2 as TimeIntegratorTrait<T>>::workspace_size(integrator),
            TimeIntegratorEnum::RungeKutta4(integrator, _) => <RungeKutta4 as TimeIntegratorTrait<T>>::workspace_size(integrator),
        }
    }
}
