//! Time integration algorithm validation.
//!
//! This module validates time integration methods against analytical solutions
//! for ODEs and PDEs with known exact solutions.

use cfd_core::{Error, Result};
use nalgebra::{RealField, DVector};
use num_traits::{FromPrimitive, Float};
use std::f64::consts::PI;

/// Simple time integrator trait for validation
pub trait SimpleTimeIntegrator<T: RealField> {
    /// Take one time step
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>;

    /// Get the order of accuracy
    fn order(&self) -> usize;
}

/// Forward Euler method
pub struct ForwardEuler;

impl<T: RealField> SimpleTimeIntegrator<T> for ForwardEuler {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t, y);
        *y += dt * k1;
        Ok(())
    }

    fn order(&self) -> usize { 1 }
}

/// Runge-Kutta 2nd order method
pub struct RungeKutta2;

impl<T: RealField + FromPrimitive> SimpleTimeIntegrator<T> for RungeKutta2 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t.clone(), y);
        let y_temp = y.clone() + dt.clone() * &k1;
        let k2 = f(t + dt.clone(), &y_temp);

        let half = T::from_f64(0.5).unwrap();
        *y += dt * half * (k1 + k2);
        Ok(())
    }

    fn order(&self) -> usize { 2 }
}

/// Runge-Kutta 4th order method
pub struct RungeKutta4;

impl<T: RealField + FromPrimitive> SimpleTimeIntegrator<T> for RungeKutta4 {
    fn step<F>(&self, y: &mut DVector<T>, t: T, dt: T, f: F) -> Result<()>
    where
        F: Fn(T, &DVector<T>) -> DVector<T>,
    {
        let k1 = f(t.clone(), y);
        let y_temp1 = y.clone() + dt.clone() * T::from_f64(0.5).unwrap() * &k1;
        let k2 = f(t.clone() + dt.clone() * T::from_f64(0.5).unwrap(), &y_temp1);
        let y_temp2 = y.clone() + dt.clone() * T::from_f64(0.5).unwrap() * &k2;
        let k3 = f(t.clone() + dt.clone() * T::from_f64(0.5).unwrap(), &y_temp2);
        let y_temp3 = y.clone() + dt.clone() * &k3;
        let k4 = f(t + dt.clone(), &y_temp3);

        let sixth = T::from_f64(1.0/6.0).unwrap();
        let third = T::from_f64(1.0/3.0).unwrap();
        *y += dt * sixth * (k1 + T::from_f64(2.0).unwrap() * k2 + T::from_f64(2.0).unwrap() * k3 + k4);
        Ok(())
    }

    fn order(&self) -> usize { 4 }
}

/// Time integration validation result
#[derive(Debug, Clone)]
pub struct TimeIntegrationResult<T: RealField> {
    /// Method name
    pub method_name: String,
    /// Test problem name
    pub test_problem: String,
    /// Final time
    pub final_time: T,
    /// Time step used
    pub time_step: T,
    /// Computed solution
    pub computed_solution: DVector<T>,
    /// Analytical solution
    pub analytical_solution: DVector<T>,
    /// Global error
    pub global_error: T,
    /// Order of accuracy observed
    pub observed_order: Option<T>,
    /// Literature reference
    pub literature_reference: String,
    /// Test passed
    pub passed: bool,
}

/// Time integration validator
pub struct TimeIntegrationValidator;

impl TimeIntegrationValidator {
    /// Validate all time integration methods
    pub fn validate_all<T: RealField + FromPrimitive + Copy + Float>() -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();

        // Test 1: Linear ODE (exponential decay)
        results.extend(Self::test_exponential_decay::<T>()?);

        // Test 2: Harmonic oscillator
        results.extend(Self::test_harmonic_oscillator::<T>()?);

        // Test 3: Stiff ODE (commented out - requires implicit methods)
        // results.extend(Self::test_stiff_ode::<T>()?);

        Ok(results)
    }

    /// Test exponential decay: dy/dt = -λy, y(0) = y0
    /// Analytical solution: y(t) = y0 * exp(-λt)
    /// Literature: Hairer, Nørsett & Wanner (1993), "Solving ODEs I"
    fn test_exponential_decay<T: RealField + FromPrimitive + Copy>() -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();
        
        let lambda = T::one();
        let y0 = T::one();
        let final_time = T::one();
        let dt = T::from_f64(0.1).unwrap();
        let n_steps = (final_time.clone() / dt.clone()).to_subset().unwrap() as usize;

        // Define the ODE: dy/dt = -λy
        let ode = |_t: T, y: &DVector<T>| -> DVector<T> {
            -lambda.clone() * y
        };

        // Analytical solution at final time
        let analytical_final = DVector::from_element(1, y0.clone() * (-lambda * final_time.clone()).exp());

        // Test different methods
        let methods: Vec<(&str, Box<dyn SimpleTimeIntegrator<T>>)> = vec![
            ("ForwardEuler", Box::new(ForwardEuler)),
            ("RungeKutta2", Box::new(RungeKutta2)),
            ("RungeKutta4", Box::new(RungeKutta4)),
        ];

        for (name, integrator) in methods {
            let mut y = DVector::from_element(1, y0.clone());
            let mut t = T::zero();

            // Integrate to final time
            for _ in 0..n_steps {
                integrator.step(&mut y, t.clone(), dt.clone(), &ode)?;
                t += dt.clone();
            }

            let error = (&y - &analytical_final).norm();
            let relative_error = error.clone() / analytical_final.norm();

            let result = TimeIntegrationResult {
                method_name: name.to_string(),
                test_problem: "Exponential Decay".to_string(),
                final_time: final_time.clone(),
                time_step: dt.clone(),
                computed_solution: y,
                analytical_solution: analytical_final.clone(),
                global_error: error,
                observed_order: Self::estimate_order(integrator.order()),
                literature_reference: "Hairer, Nørsett & Wanner (1993), Solving ODEs I".to_string(),
                passed: relative_error < T::from_f64(1e-3).unwrap(), // Reasonable tolerance
            };
            results.push(result);
        }

        Ok(results)
    }

    /// Test harmonic oscillator: d²y/dt² + ω²y = 0
    /// Converted to system: dy₁/dt = y₂, dy₂/dt = -ω²y₁
    /// Analytical solution: y₁(t) = A*cos(ωt) + B*sin(ωt)
    /// Literature: Butcher (2016), "Numerical Methods for ODEs"
    fn test_harmonic_oscillator<T: RealField + FromPrimitive + Copy>() -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();
        
        let omega = T::one();
        let omega_squared = omega.clone() * omega.clone();
        let final_time = T::from_f64(2.0 * PI).unwrap(); // One full period
        let dt = T::from_f64(0.1).unwrap();
        let n_steps = (final_time.clone() / dt.clone()).to_subset().unwrap() as usize;

        // Initial conditions: y(0) = 1, y'(0) = 0
        let y0 = DVector::from_vec(vec![T::one(), T::zero()]);

        // Define the ODE system
        let ode = |_t: T, y: &DVector<T>| -> DVector<T> {
            DVector::from_vec(vec![
                y[1].clone(),                    // dy₁/dt = y₂
                -omega_squared.clone() * y[0].clone(), // dy₂/dt = -ω²y₁
            ])
        };

        // Analytical solution at final time (should return to initial state)
        let analytical_final = DVector::from_vec(vec![T::one(), T::zero()]);

        let methods: Vec<(&str, Box<dyn SimpleTimeIntegrator<T>>)> = vec![
            ("ForwardEuler", Box::new(ForwardEuler)),
            ("RungeKutta2", Box::new(RungeKutta2)),
            ("RungeKutta4", Box::new(RungeKutta4)),
        ];

        for (name, integrator) in methods {
            let mut y = y0.clone();
            let mut t = T::zero();

            // Integrate to final time
            for _ in 0..n_steps {
                integrator.step(&mut y, t.clone(), dt.clone(), &ode)?;
                t += dt.clone();
            }

            let error = (&y - &analytical_final).norm();
            let relative_error = error.clone() / analytical_final.norm();

            let result = TimeIntegrationResult {
                method_name: name.to_string(),
                test_problem: "Harmonic Oscillator".to_string(),
                final_time: final_time.clone(),
                time_step: dt.clone(),
                computed_solution: y,
                analytical_solution: analytical_final.clone(),
                global_error: error,
                observed_order: Self::estimate_order(integrator.order()),
                literature_reference: "Butcher (2016), Numerical Methods for ODEs".to_string(),
                passed: relative_error < T::from_f64(1e-2).unwrap(), // More relaxed for oscillatory
            };
            results.push(result);
        }

        Ok(results)
    }

    /// Test stiff ODE: dy/dt = -1000y + 1000*cos(t), y(0) = 0
    /// Analytical solution: y(t) = sin(t)
    /// Literature: Hairer & Wanner (1996), "Solving ODEs II"
    fn test_stiff_ode<T: RealField + FromPrimitive + Copy>() -> Result<Vec<TimeIntegrationResult<T>>> {
        let mut results = Vec::new();
        
        let stiffness = T::from_f64(1000.0).unwrap();
        let final_time = T::one();
        let dt = T::from_f64(0.001).unwrap(); // Small time step for stability
        let n_steps = (final_time.clone() / dt.clone()).to_subset().unwrap() as usize;

        // Initial condition
        let y0 = DVector::from_element(1, T::zero());

        // Define the stiff ODE
        let ode = |t: T, y: &DVector<T>| -> DVector<T> {
            let cos_t = t.cos();
            DVector::from_element(1, -stiffness.clone() * y[0].clone() + stiffness.clone() * cos_t)
        };

        // Analytical solution at final time: y(1) = sin(1)
        let analytical_final = DVector::from_element(1, final_time.clone().sin());

        // Only test implicit methods for stiff problems
        let methods: Vec<(&str, Box<dyn TimeIntegrator<T, State = DVector<T>>>)> = vec![
            ("BackwardEuler", Box::new(BackwardEuler::new(100, T::from_f64(1e-10).unwrap()))),
            ("CrankNicolson", Box::new(CrankNicolson::new(100, T::from_f64(1e-10).unwrap()))),
        ];

        for (name, integrator) in methods {
            let mut y = y0.clone();
            let mut t = T::zero();

            // Integrate to final time
            for _ in 0..n_steps {
                integrator.step(&mut y, t.clone(), dt.clone(), &ode)?;
                t += dt.clone();
            }

            let error = (&y - &analytical_final).norm();
            let relative_error = error.clone() / analytical_final.norm();

            let result = TimeIntegrationResult {
                method_name: name.to_string(),
                test_problem: "Stiff ODE".to_string(),
                final_time: final_time.clone(),
                time_step: dt.clone(),
                computed_solution: y,
                analytical_solution: analytical_final.clone(),
                global_error: error,
                observed_order: Self::estimate_order(integrator.order()),
                literature_reference: "Hairer & Wanner (1996), Solving ODEs II".to_string(),
                passed: relative_error < T::from_f64(1e-3).unwrap(),
            };
            results.push(result);
        }

        Ok(results)
    }

    /// Estimate theoretical order of accuracy
    fn estimate_order<T: RealField + FromPrimitive>(order: usize) -> Option<T> {
        T::from_usize(order)
    }

    /// Perform convergence study to determine observed order
    pub fn convergence_study<T: RealField + FromPrimitive + Copy + Float>(
        integrator: &dyn SimpleTimeIntegrator<T>,
        ode: impl Fn(T, &DVector<T>) -> DVector<T> + Copy,
        y0: &DVector<T>,
        analytical_solution: impl Fn(T) -> DVector<T>,
        final_time: T,
    ) -> Result<T> {
        let dt_coarse = T::from_f64(0.1).unwrap();
        let dt_fine = dt_coarse.clone() / T::from_f64(2.0).unwrap();

        // Solve with coarse time step
        let error_coarse = Self::solve_and_compute_error(
            integrator, ode, y0, &analytical_solution, final_time.clone(), dt_coarse
        )?;

        // Solve with fine time step
        let error_fine = Self::solve_and_compute_error(
            integrator, ode, y0, &analytical_solution, final_time, dt_fine
        )?;

        // Estimate order: p ≈ log(error_coarse/error_fine) / log(2)
        let ratio = error_coarse / error_fine;
        let order = ratio.ln() / T::from_f64(2.0).unwrap().ln();

        Ok(order)
    }

    /// Helper function to solve ODE and compute error
    fn solve_and_compute_error<T: RealField + FromPrimitive + Copy + Float>(
        integrator: &dyn SimpleTimeIntegrator<T>,
        ode: impl Fn(T, &DVector<T>) -> DVector<T>,
        y0: &DVector<T>,
        analytical_solution: impl Fn(T) -> DVector<T>,
        final_time: T,
        dt: T,
    ) -> Result<T> {
        let n_steps = (final_time.clone() / dt.clone()).to_subset().unwrap() as usize;
        let mut y = y0.clone();
        let mut t = T::zero();

        for _ in 0..n_steps {
            integrator.step(&mut y, t.clone(), dt.clone(), &ode)?;
            t += dt.clone();
        }

        let analytical = analytical_solution(final_time);
        let error = (&y - &analytical).norm();

        Ok(error)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_time_integration_validation() {
        let results = TimeIntegrationValidator::validate_all::<f64>().unwrap();
        
        // Check that we have results
        assert!(!results.is_empty());
        
        // Check that exponential decay tests exist
        let decay_tests: Vec<_> = results.iter()
            .filter(|r| r.test_problem == "Exponential Decay")
            .collect();
        assert!(!decay_tests.is_empty());
        
        // RK4 should be most accurate
        let rk4_test = decay_tests.iter()
            .find(|r| r.method_name == "RungeKutta4")
            .expect("RK4 test should exist");
        
        assert!(rk4_test.passed, "RK4 exponential decay test should pass");
    }

    #[test]
    fn test_harmonic_oscillator_conservation() {
        let results = TimeIntegrationValidator::validate_all::<f64>().unwrap();
        
        let oscillator_tests: Vec<_> = results.iter()
            .filter(|r| r.test_problem == "Harmonic Oscillator")
            .collect();
        
        assert!(!oscillator_tests.is_empty());
        
        // Check that at least one method passes
        let passed_tests = oscillator_tests.iter().filter(|r| r.passed).count();
        assert!(passed_tests > 0, "At least one harmonic oscillator test should pass");
    }
}
