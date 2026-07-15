//! Traits for time-stepping methods.

use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};

/// Leto-backed ODE state vector used by time-stepping methods.
pub type TimeState<T> = Array1<T>;

/// Leto-backed dense matrix used by time-stepping methods.
pub type TimeMatrix<T> = Array2<T>;

pub(crate) fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

pub(crate) fn zero<T: NumericElement>() -> T {
    <T as NumericElement>::ZERO
}

pub(crate) fn one<T: NumericElement>() -> T {
    <T as NumericElement>::ONE
}

#[cfg(test)]
pub(crate) fn state_from_vec<T>(values: Vec<T>) -> TimeState<T> {
    TimeState::from_shape_vec([values.len()], values)
        .expect("invariant: state length matches Leto rank-1 shape")
}

pub(crate) fn state_zeros<T: NumericElement + Clone>(len: usize) -> TimeState<T> {
    TimeState::from_elem([len], zero())
}

pub(crate) fn state_len<T>(state: &TimeState<T>) -> usize {
    state.shape()[0]
}

pub(crate) fn ensure_same_len<T>(
    lhs: &TimeState<T>,
    rhs: &TimeState<T>,
    context: &str,
) -> Result<usize> {
    let lhs_len = state_len(lhs);
    let rhs_len = state_len(rhs);
    if lhs_len != rhs_len {
        return Err(Error::InvalidInput(format!(
            "{context}: lhs state has {lhs_len} entries, rhs state has {rhs_len}"
        )));
    }
    Ok(lhs_len)
}

pub(crate) fn state_norm<T: RealField + Copy>(state: &TimeState<T>) -> T {
    let mut sum = zero();
    for i in 0..state_len(state) {
        sum += state[i] * state[i];
    }
    <T as NumericElement>::sqrt(sum)
}

pub(crate) fn add_scaled_in_place<T: RealField + Copy>(
    target: &mut TimeState<T>,
    increment: &TimeState<T>,
    scale: T,
    context: &str,
) -> Result<()> {
    let len = ensure_same_len(target, increment, context)?;
    for i in 0..len {
        target[i] += scale * increment[i];
    }
    Ok(())
}

pub(crate) fn assign_base_plus_scaled<T: RealField + Copy>(
    target: &mut TimeState<T>,
    base: &TimeState<T>,
    increment: &TimeState<T>,
    scale: T,
    context: &str,
) -> Result<()> {
    let len = ensure_same_len(base, increment, context)?;
    if state_len(target) != len {
        return Err(Error::InvalidInput(format!(
            "{context}: target state has {} entries, base state has {len}",
            state_len(target)
        )));
    }
    for i in 0..len {
        target[i] = base[i] + scale * increment[i];
    }
    Ok(())
}

#[cfg(test)]
pub(crate) fn state_scale<T: RealField + Copy>(state: &TimeState<T>, scale: T) -> TimeState<T> {
    state_from_vec((0..state_len(state)).map(|i| state[i] * scale).collect())
}

#[cfg(test)]
pub(crate) fn state_neg<T: RealField + Copy>(state: &TimeState<T>) -> TimeState<T> {
    state_scale(state, -one::<T>())
}

/// Time-stepping method for ODE integration du/dt = f(t,u)
pub trait TimeStepper<T: RealField + Copy> {
    /// Advance solution by one time step
    ///
    /// # Arguments
    /// * `f` - Right-hand side function f(t, u)
    /// * `t` - Current time
    /// * `u` - Current solution vector
    /// * `dt` - Time step size
    ///
    /// # Returns
    /// New solution vector at t + dt
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>;

    /// Get the order of accuracy of this method
    fn order(&self) -> usize;

    /// Get the number of stages (function evaluations per step)
    fn stages(&self) -> usize;

    /// Check if this is an explicit method
    fn is_explicit(&self) -> bool {
        true // Default to explicit
    }

    /// Get method stability region information
    fn stability_region(&self) -> Option<&str> {
        None
    }
}

/// Controller for adaptive time stepping
pub trait TimeStepController<T: RealField + Copy> {
    /// Compute optimal time step based on error estimate
    ///
    /// # Arguments
    /// * `error_estimate` - Estimated local truncation error
    /// * `current_dt` - Current time step size
    /// * `tolerance` - Desired error tolerance
    ///
    /// # Returns
    /// Recommended new time step size
    fn adapt_step(&self, error_estimate: T, current_dt: T, tolerance: T) -> T;

    /// Check if solution should be accepted
    ///
    /// # Arguments
    /// * `error_estimate` - Estimated local truncation error
    /// * `tolerance` - Desired error tolerance
    ///
    /// # Returns
    /// true if step should be accepted, false if it should be rejected
    fn accept_step(&self, error_estimate: T, tolerance: T) -> bool {
        error_estimate <= tolerance
    }
}

/// Embedded Runge-Kutta method providing error estimation
pub trait EmbeddedMethod<T: RealField + Copy>: TimeStepper<T> {
    /// Take a step and return both solution and error estimate
    ///
    /// # Arguments
    /// * `f` - Right-hand side function f(t, u)
    /// * `t` - Current time
    /// * `u` - Current solution vector
    /// * `dt` - Time step size
    ///
    /// # Returns
    /// Tuple of (solution, error_estimate)
    fn embedded_step<F>(
        &self,
        f: F,
        t: T,
        u: &TimeState<T>,
        dt: T,
    ) -> Result<(TimeState<T>, TimeState<T>)>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>;
}
