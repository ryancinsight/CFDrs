//! Time integration methods for validation.

use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};
use leto::Array1;
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

/// Leto-backed ODE state vector used by validation time integrators.
pub type State<T> = Array1<T>;

/// Create a Leto rank-1 state from a contiguous vector.
pub(crate) fn state_from_vec<T>(values: Vec<T>) -> State<T> {
    State::from_shape_vec([values.len()], values)
        .expect("invariant: state length matches Leto rank-1 shape")
}

/// Create a Leto rank-1 state filled with one value.
pub(crate) fn state_from_elem<T: Copy>(len: usize, value: T) -> State<T> {
    State::from_elem([len], value)
}

/// Create a zero-filled Leto rank-1 state.
pub(crate) fn state_zeros<T: RealField + Copy>(len: usize) -> State<T> {
    State::from_elem([len], scalar::zero())
}

/// Return the state length.
pub(crate) fn state_len<T>(state: &State<T>) -> usize {
    state.shape()[0]
}

fn ensure_same_len<T>(lhs: &State<T>, rhs: &State<T>) -> Result<usize> {
    let lhs_len = state_len(lhs);
    let rhs_len = state_len(rhs);
    if lhs_len != rhs_len {
        return Err(Error::InvalidInput(format!(
            "Time-integration state dimension mismatch: lhs has {lhs_len}, rhs has {rhs_len}"
        )));
    }
    Ok(lhs_len)
}

fn add_scaled_in_place<T: RealField + Copy>(
    target: &mut State<T>,
    increment: &State<T>,
    scale: T,
) -> Result<()> {
    let len = ensure_same_len(target, increment)?;
    for i in 0..len {
        target[i] += increment[i] * scale;
    }
    Ok(())
}

fn assign_base_plus_scaled<T: RealField + Copy>(
    target: &mut State<T>,
    base: &State<T>,
    increment: &State<T>,
    scale: T,
) -> Result<()> {
    let len = ensure_same_len(base, increment)?;
    if state_len(target) != len {
        return Err(Error::InvalidInput(format!(
            "Time-integration workspace dimension mismatch: workspace has {}, state has {len}",
            state_len(target)
        )));
    }
    for i in 0..len {
        target[i] = base[i] + increment[i] * scale;
    }
    Ok(())
}

/// Time integrator trait for validation
pub trait TimeIntegratorTrait<T: RealField + Copy> {
    /// Take one time step with workspace buffer for zero-copy operation
    fn step<F>(&self, y: &mut State<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>;

    /// Take one time step with pre-allocated workspace for zero-copy performance
    fn step_with_workspace<F>(
        &self,
        y: &mut State<T>,
        workspace: &mut [State<T>],
        t: T,
        dt: T,
        f: F,
    ) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>;

    /// Get the order of accuracy
    fn order(&self) -> usize;

    /// Get required workspace size for zero-copy operation
    fn workspace_size(&self) -> usize;
}

/// Forward Euler time integrator (first-order explicit method)
pub struct ForwardEuler;

impl<T: RealField + FloatElement + Copy> TimeIntegratorTrait<T> for ForwardEuler {
    fn step<F>(&self, y: &mut State<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        let k1 = f(t, y);
        add_scaled_in_place(y, &k1, dt)
    }

    fn step_with_workspace<F>(
        &self,
        y: &mut State<T>,
        _workspace: &mut [State<T>],
        t: T,
        dt: T,
        f: F,
    ) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
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

impl<T: RealField + FloatElement + Copy> TimeIntegratorTrait<T> for RungeKutta2 {
    fn step<F>(&self, y: &mut State<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        // Default implementation using internal workspace allocation
        let mut workspace =
            vec![state_zeros(state_len(y)); <Self as TimeIntegratorTrait<T>>::workspace_size(self)];
        self.step_with_workspace(y, &mut workspace, t, dt, f)
    }

    fn step_with_workspace<F>(
        &self,
        y: &mut State<T>,
        workspace: &mut [State<T>],
        t: T,
        dt: T,
        f: F,
    ) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        if workspace.len() < <Self as TimeIntegratorTrait<T>>::workspace_size(self) {
            return Err(cfd_core::error::Error::InvalidInput(
                "Insufficient workspace size".to_string(),
            ));
        }

        let k1 = f(t, y);

        // True zero-copy: Use pre-allocated workspace buffer
        let y_intermediate = &mut workspace[0];
        assign_base_plus_scaled(y_intermediate, y, &k1, dt)?;

        let k2 = f(t + dt, y_intermediate);

        let half = scalar::from_f64::<T>(HALF);
        let scale = dt * half;
        let len = ensure_same_len(&k1, &k2)?;
        ensure_same_len(y, &k1)?;
        for i in 0..len {
            y[i] += (k1[i] + k2[i]) * scale;
        }
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

impl<T: RealField + FloatElement + Copy> TimeIntegratorTrait<T> for RungeKutta4 {
    fn step<F>(&self, y: &mut State<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        // Default implementation using internal workspace allocation
        let mut workspace =
            vec![state_zeros(state_len(y)); <Self as TimeIntegratorTrait<T>>::workspace_size(self)];
        self.step_with_workspace(y, &mut workspace, t, dt, f)
    }

    fn step_with_workspace<F>(
        &self,
        y: &mut State<T>,
        workspace: &mut [State<T>],
        t: T,
        dt: T,
        f: F,
    ) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        if workspace.len() < <Self as TimeIntegratorTrait<T>>::workspace_size(self) {
            return Err(cfd_core::error::Error::InvalidInput(
                "Insufficient workspace size".to_string(),
            ));
        }

        let half = scalar::from_f64::<T>(HALF);

        let k1 = f(t, y);

        // True zero-copy: Use pre-allocated workspace buffers
        let y_k1 = &mut workspace[0];
        assign_base_plus_scaled(y_k1, y, &k1, dt * half)?;

        let k2 = f(t + dt * half, y_k1);

        let y_k2 = &mut workspace[1];
        assign_base_plus_scaled(y_k2, y, &k2, dt * half)?;

        let k3 = f(t + dt * half, y_k2);

        let y_k3 = &mut workspace[2];
        assign_base_plus_scaled(y_k3, y, &k3, dt)?;

        let k4 = f(t + dt, y_k3);

        let sixth = scalar::from_f64::<T>(ONE_SIXTH);
        let two = scalar::from_f64::<T>(TWO);
        let scale = dt * sixth;
        let len = ensure_same_len(&k1, &k2)?;
        ensure_same_len(&k1, &k3)?;
        ensure_same_len(&k1, &k4)?;
        ensure_same_len(y, &k1)?;
        for i in 0..len {
            y[i] += (k1[i] + k2[i] * two + k3[i] * two + k4[i]) * scale;
        }
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

impl<T: RealField + FloatElement + Copy> TimeIntegratorTrait<T> for TimeIntegratorEnum<T> {
    fn step<F>(&self, y: &mut State<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        match self {
            TimeIntegratorEnum::ForwardEuler(integrator, _) => integrator.step(y, t, dt, f),
            TimeIntegratorEnum::RungeKutta2(integrator, _) => integrator.step(y, t, dt, f),
            TimeIntegratorEnum::RungeKutta4(integrator, _) => integrator.step(y, t, dt, f),
        }
    }

    fn step_with_workspace<F>(
        &self,
        y: &mut State<T>,
        workspace: &mut [State<T>],
        t: T,
        dt: T,
        f: F,
    ) -> Result<()>
    where
        F: Fn(T, &State<T>) -> State<T>,
    {
        match self {
            TimeIntegratorEnum::ForwardEuler(integrator, _) => {
                integrator.step_with_workspace(y, workspace, t, dt, f)
            }
            TimeIntegratorEnum::RungeKutta2(integrator, _) => {
                integrator.step_with_workspace(y, workspace, t, dt, f)
            }
            TimeIntegratorEnum::RungeKutta4(integrator, _) => {
                integrator.step_with_workspace(y, workspace, t, dt, f)
            }
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
            TimeIntegratorEnum::ForwardEuler(integrator, _) => {
                <ForwardEuler as TimeIntegratorTrait<T>>::workspace_size(integrator)
            }
            TimeIntegratorEnum::RungeKutta2(integrator, _) => {
                <RungeKutta2 as TimeIntegratorTrait<T>>::workspace_size(integrator)
            }
            TimeIntegratorEnum::RungeKutta4(integrator, _) => {
                <RungeKutta4 as TimeIntegratorTrait<T>>::workspace_size(integrator)
            }
        }
    }
}
