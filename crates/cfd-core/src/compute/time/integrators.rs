//! Time integration schemes.

use crate::error::Result;
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

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
    type State = Array1<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let derivative = f(t, state);
        add_scaled_in_place(state, &derivative, dt)?;
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

impl<T: RealField + FloatElement + Copy> TimeIntegrator<T> for RungeKutta2 {
    type State = Array1<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let half = <T as FloatElement>::from_f64(0.5);

        let k1 = f(t, state);
        let mut state_buffer = state.clone();
        add_scaled_in_place(&mut state_buffer, &k1, dt * half)?;

        let k2 = f(t + dt * half, &state_buffer);
        add_scaled_in_place(state, &k2, dt)?;

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

impl<T: RealField + FloatElement + Copy> TimeIntegrator<T> for RungeKutta4 {
    type State = Array1<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        let two = <T as FloatElement>::from_f64(2.0);
        let six = <T as FloatElement>::from_f64(6.0);
        let half = <T as FloatElement>::from_f64(0.5);

        let k1 = f(t, state);

        let mut state_buffer = state.clone();
        add_scaled_in_place(&mut state_buffer, &k1, dt * half)?;
        let k2 = f(t + dt * half, &state_buffer);

        state_buffer.clone_from(state);
        add_scaled_in_place(&mut state_buffer, &k2, dt * half)?;
        let k3 = f(t + dt * half, &state_buffer);

        state_buffer.clone_from(state);
        add_scaled_in_place(&mut state_buffer, &k3, dt)?;
        let k4 = f(t + dt, &state_buffer);

        // y_{n+1} = y_n + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        add_scaled_in_place(state, &k1, dt / six)?;
        add_scaled_in_place(state, &k2, dt * two / six)?;
        add_scaled_in_place(state, &k3, dt * two / six)?;
        add_scaled_in_place(state, &k4, dt / six)?;

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

impl<T: RealField + FloatElement + Copy> Default for BackwardEuler<T> {
    fn default() -> Self {
        Self {
            tolerance: <T as FloatElement>::from_f64(1e-10),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FloatElement + Copy> BackwardEuler<T> {
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
    type State = Array1<T>;

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
            let next_state = add_scaled(&previous_state, &f_val, dt)?;

            // Check convergence
            let error = difference_norm(&next_state, &next_candidate)?;
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

impl<T: RealField + FloatElement + Copy> Default for CrankNicolson<T> {
    fn default() -> Self {
        Self {
            tolerance: <T as FloatElement>::from_f64(1e-10),
            max_iterations: 100,
        }
    }
}

impl<T: RealField + FloatElement + Copy> CrankNicolson<T> {
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

impl<T: RealField + FloatElement + Copy> TimeIntegrator<T> for CrankNicolson<T> {
    type State = Array1<T>;

    fn step<F>(&self, state: &mut Self::State, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &Self::State) -> Self::State,
    {
        // Implement Crank-Nicolson using fixed-point iteration
        // Solve: y_{n+1} = y_n + dt/2 * (f(t_n, y_n) + f(t_{n+1}, y_{n+1}))

        let mut next_candidate = state.clone();
        let previous_state = state.clone();
        let t_next = t + dt;
        let half = <T as FloatElement>::from_f64(0.5);
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
            let f_average = add(&f_previous, &f_current)?;
            let next_state = add_scaled(&previous_state, &f_average, half_dt)?;

            // Check convergence
            let error = difference_norm(&next_state, &next_candidate)?;
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

fn ensure_same_len<T>(left: &Array1<T>, right: &Array1<T>, context: &str) -> Result<()> {
    if left.size() == right.size() {
        Ok(())
    } else {
        Err(crate::error::Error::InvalidConfiguration(format!(
            "{context}: expected equal state lengths, got {} and {}",
            left.size(),
            right.size()
        )))
    }
}

fn contiguous_slice<T>(array: &Array1<T>) -> &[T] {
    array
        .as_slice()
        .expect("invariant: cfd-core time states use contiguous Leto Array1 storage")
}

fn contiguous_slice_mut<T>(array: &mut Array1<T>) -> &mut [T] {
    array
        .as_slice_mut()
        .expect("invariant: cfd-core time states use contiguous Leto Array1 storage")
}

fn add_scaled_in_place<T>(state: &mut Array1<T>, derivative: &Array1<T>, scale: T) -> Result<()>
where
    T: NumericElement,
{
    ensure_same_len(state, derivative, "time integration scaled update")?;

    for (slot, &increment) in contiguous_slice_mut(state)
        .iter_mut()
        .zip(contiguous_slice(derivative))
    {
        *slot += increment * scale;
    }
    Ok(())
}

fn add_scaled<T>(state: &Array1<T>, derivative: &Array1<T>, scale: T) -> Result<Array1<T>>
where
    T: NumericElement,
{
    ensure_same_len(state, derivative, "time integration scaled sum")?;
    let values = contiguous_slice(state)
        .iter()
        .zip(contiguous_slice(derivative))
        .map(|(&value, &increment)| value + increment * scale)
        .collect();
    Array1::from_vec([state.size()], values).map_err(|error| {
        crate::error::Error::InvalidConfiguration(format!(
            "time integration scaled sum shape: {error}"
        ))
    })
}

fn add<T>(left: &Array1<T>, right: &Array1<T>) -> Result<Array1<T>>
where
    T: NumericElement,
{
    ensure_same_len(left, right, "time integration state sum")?;
    let values = contiguous_slice(left)
        .iter()
        .zip(contiguous_slice(right))
        .map(|(&lhs, &rhs)| lhs + rhs)
        .collect();
    Array1::from_vec([left.size()], values).map_err(|error| {
        crate::error::Error::InvalidConfiguration(format!(
            "time integration state sum shape: {error}"
        ))
    })
}

fn difference_norm<T>(left: &Array1<T>, right: &Array1<T>) -> Result<T>
where
    T: NumericElement,
{
    ensure_same_len(left, right, "time integration convergence check")?;
    let squared = contiguous_slice(left)
        .iter()
        .zip(contiguous_slice(right))
        .fold(<T as NumericElement>::ZERO, |sum, (&lhs, &rhs)| {
            let diff = lhs - rhs;
            sum + diff * diff
        });
    Ok(squared.sqrt())
}

#[cfg(test)]
mod tests {
    use super::{
        BackwardEuler, CrankNicolson, ForwardEuler, RungeKutta2, RungeKutta4, TimeIntegrator,
    };
    use crate::error::{ConvergenceErrorKind, Error};
    use leto::Array1;

    fn state_from(values: Vec<f64>) -> Array1<f64> {
        Array1::from_vec([values.len()], values).expect("test state shape matches supplied values")
    }

    #[test]
    fn forward_euler_applies_constant_derivative() {
        let mut state = state_from(vec![1.0, -2.0]);

        ForwardEuler
            .step(&mut state, 0.0, 0.25, |_, _| state_from(vec![2.0, -4.0]))
            .expect("constant derivative step is well-defined");

        assert_eq!(
            state
                .as_slice()
                .expect("test state uses contiguous Leto Array1 storage"),
            &[1.5, -3.0]
        );
        assert_eq!(
            <ForwardEuler as TimeIntegrator<f64>>::order(&ForwardEuler),
            1
        );
        assert!(<ForwardEuler as TimeIntegrator<f64>>::is_explicit(
            &ForwardEuler
        ));
    }

    #[test]
    fn runge_kutta2_applies_midpoint_update() {
        let mut state = state_from(vec![1.0]);

        RungeKutta2
            .step(&mut state, 0.0, 0.5, |t, y| state_from(vec![t + y[[0]]]))
            .expect("linear midpoint update is well-defined");

        assert_eq!(state[[0]], 1.75);
        assert_eq!(<RungeKutta2 as TimeIntegrator<f64>>::order(&RungeKutta2), 2);
        assert!(<RungeKutta2 as TimeIntegrator<f64>>::is_explicit(
            &RungeKutta2
        ));
    }

    #[test]
    fn runge_kutta4_matches_fourth_order_taylor_for_exponential_growth() {
        let mut state = state_from(vec![1.0]);

        RungeKutta4
            .step(&mut state, 0.0, 0.1, |_, y| y.clone())
            .expect("autonomous linear step is well-defined");

        let expected: f64 = 1.0 + 0.1 + 0.005 + (0.001 / 6.0) + (0.0001 / 24.0);
        let tolerance: f64 = 16.0 * f64::EPSILON;
        assert!(
            (state[[0]] - expected).abs() <= tolerance,
            "state = {}, expected = {}, tolerance = {}",
            state[[0]],
            expected,
            tolerance
        );
        assert_eq!(<RungeKutta4 as TimeIntegrator<f64>>::order(&RungeKutta4), 4);
        assert!(<RungeKutta4 as TimeIntegrator<f64>>::is_explicit(
            &RungeKutta4
        ));
    }

    #[test]
    fn implicit_integrator_defaults_are_value_semantic() {
        let backward = BackwardEuler::<f64>::default();
        assert_eq!(backward.tolerance, 1e-10);
        assert_eq!(backward.max_iterations, 100);
        assert_eq!(
            <BackwardEuler<f64> as TimeIntegrator<f64>>::order(&backward),
            1
        );
        assert!(!<BackwardEuler<f64> as TimeIntegrator<f64>>::is_explicit(
            &backward
        ));

        let crank = CrankNicolson::<f64>::default();
        assert_eq!(crank.tolerance, 1e-10);
        assert_eq!(crank.max_iterations, 100);
        assert_eq!(
            <CrankNicolson<f64> as TimeIntegrator<f64>>::order(&crank),
            2
        );
        assert!(!<CrankNicolson<f64> as TimeIntegrator<f64>>::is_explicit(
            &crank
        ));
    }

    #[test]
    fn implicit_integrators_reject_zero_iterations() {
        let mut state = state_from(vec![1.0]);
        let backward = BackwardEuler::<f64>::new().with_max_iterations(0);
        let error = backward
            .step(&mut state, 0.0, 0.1, |_, y| y.clone())
            .expect_err("zero iterations must be rejected");
        assert!(
            matches!(error, Error::InvalidConfiguration(message) if message.contains("BackwardEuler"))
        );

        let mut state = state_from(vec![1.0]);
        let crank = CrankNicolson::<f64>::new().with_max_iterations(0);
        let error = crank
            .step(&mut state, 0.0, 0.1, |_, y| y.clone())
            .expect_err("zero iterations must be rejected");
        assert!(
            matches!(error, Error::InvalidConfiguration(message) if message.contains("CrankNicolson"))
        );
    }

    #[test]
    fn backward_euler_reports_nonconvergence_value_semantically() {
        let mut state = state_from(vec![1.0]);
        let integrator = BackwardEuler::<f64>::new()
            .with_tolerance(1e-14)
            .with_max_iterations(1);

        let error = integrator
            .step(&mut state, 0.0, 0.1, |_, y| y.clone())
            .expect_err("one iteration cannot satisfy this implicit equation");

        assert!(matches!(
            error,
            Error::Convergence(ConvergenceErrorKind::MaxIterationsExceeded { max: 1 })
        ));
    }
}
