//! Time integration validation tests.

use super::integrators::{ForwardEuler, RungeKutta2, RungeKutta4, TimeIntegratorTrait};
use super::results::TimeIntegrationResult;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::{FromPrimitive, ToPrimitive};
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
    pub fn validate_all<T: RealField + Copy + FromPrimitive + ToPrimitive>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        // Test exponential decay
        results.extend(Self::validate_exponential_decay()?);

        // Test harmonic oscillator
        results.extend(Self::validate_harmonic_oscillator()?);

        Ok(results)
    }

    /// Validate exponential decay: dy/dt = -λy
    fn validate_exponential_decay<T: RealField + Copy + FromPrimitive + ToPrimitive>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        let lambda = T::from_f64(DECAY_LAMBDA).unwrap_or_else(T::zero);
        let dt = T::from_f64(TIME_STEP_VALIDATION).unwrap_or_else(T::zero);
        let final_time = T::one();
        let n_steps = (final_time / dt).to_subset().unwrap_or(100.0) as usize;

        // Initial condition
        let y0 = DVector::from_element(1, T::one());

        // ODE function
        let f = |_t: T, y: &DVector<T>| -> DVector<T> { y * (-lambda) };

        // Test Forward Euler
        let mut y = y0.clone();
        let euler = ForwardEuler;
        for _ in 0..n_steps {
            euler.step(&mut y, T::zero(), dt, &f)?;
        }
        let exact = (-lambda * final_time).exp();
        let error = (y[0] - exact).abs();

        results.push(TimeIntegrationResult::create(
            "ForwardEuler".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < T::from_f64(ERROR_THRESHOLD).unwrap_or_else(T::zero),
            1,
        ));

        // Test RK2
        let mut y = y0.clone();
        let rk2 = RungeKutta2;
        for _ in 0..n_steps {
            rk2.step(&mut y, T::zero(), dt, &f)?;
        }
        let error = (y[0] - exact).abs();

        results.push(TimeIntegrationResult::create(
            "RungeKutta2".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < T::from_f64(ERROR_THRESHOLD).unwrap_or_else(T::zero),
            2,
        ));

        // Test RK4
        let mut y = y0.clone();
        let rk4 = RungeKutta4;
        for _ in 0..n_steps {
            rk4.step(&mut y, T::zero(), dt, &f)?;
        }
        let error = (y[0] - exact).abs();

        results.push(TimeIntegrationResult::create(
            "RungeKutta4".to_string(),
            "ExponentialDecay".to_string(),
            dt,
            final_time,
            n_steps,
            error,
            error < T::from_f64(ERROR_THRESHOLD).unwrap_or_else(T::zero),
            4,
        ));

        Ok(results)
    }

    /// Validate harmonic oscillator: d²y/dt² + ω²y = 0
    fn validate_harmonic_oscillator<T: RealField + Copy + FromPrimitive + ToPrimitive>(
    ) -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        let omega = T::from_f64(OSCILLATOR_OMEGA).unwrap_or_else(T::zero);
        let dt = T::from_f64(TIME_STEP_VALIDATION).unwrap_or_else(T::zero);
        let final_time = T::from_f64(2.0 * PI).unwrap_or_else(T::zero); // One period
        let n_steps = (final_time / dt).to_subset().unwrap_or(628.0) as usize;

        // Initial conditions: y(0) = 1, y'(0) = 0
        let y0 = DVector::from_vec(vec![T::one(), T::zero()]);

        // Convert to first-order system
        let f = |_t: T, y: &DVector<T>| -> DVector<T> {
            DVector::from_vec(vec![y[1], -omega * omega * y[0]])
        };

        // Test each integrator
        let integrators: Vec<(&str, Box<dyn Fn(&mut DVector<T>, T, T) -> Result<()>>)> = vec![
            (
                "ForwardEuler",
                Box::new(|y: &mut DVector<T>, t, dt| ForwardEuler.step(y, t, dt, &f)),
            ),
            (
                "RungeKutta2",
                Box::new(|y: &mut DVector<T>, t, dt| RungeKutta2.step(y, t, dt, &f)),
            ),
            (
                "RungeKutta4",
                Box::new(|y: &mut DVector<T>, t, dt| RungeKutta4.step(y, t, dt, &f)),
            ),
        ];

        for (name, integrator_fn) in integrators {
            let mut y = y0.clone();
            let mut t = T::zero();

            for _ in 0..n_steps {
                integrator_fn(&mut y, t, dt)?;
                t = t + dt;
            }

            // Exact solution: y(t) = cos(ωt)
            let exact = (omega * final_time).cos();
            let error = (y[0] - exact).abs();

            results.push(TimeIntegrationResult::create(
                name.to_string(),
                "HarmonicOscillator".to_string(),
                dt,
                final_time,
                n_steps,
                error,
                error < T::from_f64(ERROR_THRESHOLD).unwrap_or_else(T::zero),
                match name {
                    "ForwardEuler" => 1,
                    "RungeKutta2" => 2,
                    "RungeKutta4" => 4,
                    _ => 0,
                },
            ));
        }

        Ok(results)
    }
}
