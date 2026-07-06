//! Time integration validation tests.

use super::integrators::{
    state_from_elem, state_from_vec, state_len, ForwardEuler, RungeKutta2, RungeKutta4, State,
    TimeIntegratorTrait,
};
use super::results::TimeIntegrationResult;
use crate::scalar;
use cfd_core::error::Result;
use eunomia::{FloatElement, RealField};
use std::f64::consts::PI;

// Validation constants
const TIME_STEP_VALIDATION: f64 = 0.01;
const DECAY_LAMBDA: f64 = 1.0;
const OSCILLATOR_OMEGA: f64 = 1.0;
const ERROR_THRESHOLD: f64 = 1e-3;

/// Time integration validator
pub struct TimeIntegrationValidator;

impl TimeIntegrationValidator {
    /// Run all validation tests
    pub fn validate_all<T: RealField + Copy + FloatElement>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        // Test exponential decay
        results.extend(Self::validate_exponential_decay()?);

        // Test harmonic oscillator
        results.extend(Self::validate_harmonic_oscillator()?);

        Ok(results)
    }

    /// Validate exponential decay: dy/dt = -λy
    fn validate_exponential_decay<T: RealField + Copy + FloatElement>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        let lambda = scalar::from_f64::<T>(DECAY_LAMBDA);
        let dt = scalar::from_f64::<T>(TIME_STEP_VALIDATION);
        let final_time = scalar::one();

        // Safe conversion with bounds checking: clamp to reasonable range [1, 1M]
        let n_steps_f64 = scalar::to_f64(final_time / dt);
        #[allow(clippy::cast_possible_truncation)] // Clamped to [1, 1M], rounded - safe for usize
        let n_steps = n_steps_f64.clamp(1.0, 1_000_000.0).round() as usize;

        // Initial condition
        let y0 = state_from_elem(1, scalar::one());

        // ODE function
        let f = |_t: T, y: &State<T>| -> State<T> { Self::scale_state(y, -lambda) };

        // Test Forward Euler
        let mut y = y0.clone();
        let euler = ForwardEuler;
        for _ in 0..n_steps {
            euler.step(&mut y, scalar::zero(), dt, f)?;
        }
        let exact = scalar::exp(-lambda * final_time);
        let error = scalar::abs(y[0] - exact);

        results.push(TimeIntegrationResult::create(
            "ForwardEuler".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < scalar::from_f64::<T>(ERROR_THRESHOLD),
            1,
        ));

        // Test RK2
        let mut y = y0.clone();
        let rk2 = RungeKutta2;
        for _ in 0..n_steps {
            rk2.step(&mut y, scalar::zero(), dt, f)?;
        }
        let error = scalar::abs(y[0] - exact);

        results.push(TimeIntegrationResult::create(
            "RungeKutta2".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < scalar::from_f64::<T>(ERROR_THRESHOLD),
            2,
        ));

        // Test RK4
        let mut y = y0.clone();
        let rk4 = RungeKutta4;
        for _ in 0..n_steps {
            rk4.step(&mut y, scalar::zero(), dt, f)?;
        }
        let error = scalar::abs(y[0] - exact);

        results.push(TimeIntegrationResult::create(
            "RungeKutta4".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < scalar::from_f64::<T>(ERROR_THRESHOLD),
            4,
        ));

        Ok(results)
    }

    /// Validate harmonic oscillator: d²y/dt² + ω²y = 0
    fn validate_harmonic_oscillator<T: RealField + Copy + FloatElement>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        let omega = scalar::from_f64::<T>(OSCILLATOR_OMEGA);
        let dt = scalar::from_f64::<T>(TIME_STEP_VALIDATION);
        let final_time = scalar::from_f64::<T>(2.0 * PI); // One period

        // Safe conversion with bounds checking: clamp to reasonable range [1, 1M]
        let n_steps_f64 = scalar::to_f64(final_time / dt);
        #[allow(clippy::cast_possible_truncation)] // Clamped to [1, 1M], rounded - safe for usize
        let n_steps = n_steps_f64.clamp(1.0, 1_000_000.0).round() as usize;

        // Initial conditions: y(0) = 1, y'(0) = 0
        let y0 = state_from_vec(vec![scalar::one(), scalar::zero()]);

        // Convert to first-order system
        let f =
            |_t: T, y: &State<T>| -> State<T> { state_from_vec(vec![y[1], -omega * omega * y[0]]) };

        results.push(Self::validate_oscillator_integrator(
            "ForwardEuler",
            1,
            ForwardEuler,
            &y0,
            final_time,
            dt,
            n_steps,
            omega,
            f,
        )?);
        results.push(Self::validate_oscillator_integrator(
            "RungeKutta2",
            2,
            RungeKutta2,
            &y0,
            final_time,
            dt,
            n_steps,
            omega,
            f,
        )?);
        results.push(Self::validate_oscillator_integrator(
            "RungeKutta4",
            4,
            RungeKutta4,
            &y0,
            final_time,
            dt,
            n_steps,
            omega,
            f,
        )?);

        Ok(results)
    }

    fn validate_oscillator_integrator<T, I, F>(
        name: &str,
        order: usize,
        integrator: I,
        y0: &State<T>,
        final_time: T,
        dt: T,
        n_steps: usize,
        omega: T,
        f: F,
    ) -> Result<TimeIntegrationResult<T>>
    where
        T: RealField + Copy + FloatElement,
        I: TimeIntegratorTrait<T>,
        F: Copy + Fn(T, &State<T>) -> State<T>,
    {
        let mut y = y0.clone();
        let mut t = scalar::zero();

        for _ in 0..n_steps {
            integrator.step(&mut y, t, dt, f)?;
            t += dt;
        }

        // Exact solution: y(t) = cos(ωt)
        let exact = scalar::cos(omega * final_time);
        let error = scalar::abs(y[0] - exact);

        Ok(TimeIntegrationResult::create(
            name.to_string(),
            "HarmonicOscillator".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < scalar::from_f64::<T>(ERROR_THRESHOLD),
            order,
        ))
    }

    fn scale_state<T: RealField + Copy>(state: &State<T>, scale: T) -> State<T> {
        let mut values = Vec::with_capacity(state_len(state));
        for i in 0..state_len(state) {
            values.push(state[i] * scale);
        }
        state_from_vec(values)
    }
}
