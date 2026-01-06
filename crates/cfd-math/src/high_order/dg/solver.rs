//! Time integration solvers for Discontinuous Galerkin methods.
//!
//! This module provides various time integration schemes for evolving DG solutions in time,
//! including explicit, implicit, and IMEX (Implicit-Explicit) methods.

use super::{Result, DGError, DGOperator, legendre_poly};
use nalgebra::{DVector, DMatrix};
use std::time::Instant;

/// Type for the right-hand side function
pub type RhsFn<'a> = dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>> + 'a;

/// Type for the Jacobian function
pub type JacobianFn<'a> = dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>> + 'a;

/// Type of time integration method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeIntegration {
    /// Forward Euler (1st order explicit)
    ForwardEuler,
    /// Runge-Kutta 2nd order (Heun's method)
    RK2,
    /// Runge-Kutta 3rd order
    RK3,
    /// Classic 4th order Runge-Kutta
    RK4,
    /// Strong Stability Preserving RK3 (3rd order)
    SSPRK3,
    /// Adams-Bashforth 2nd order
    AB2,
    /// Adams-Bashforth 3rd order
    AB3,
    /// Implicit Euler (1st order)
    ImplicitEuler,
    /// Crank-Nicolson (2nd order)
    CrankNicolson,
    /// DIRK2 (2nd order Diagonally Implicit RK)
    DIRK2,
    /// IMEX RK2 (2nd order Implicit-Explicit)
    IMEXRK2,
    /// IMEX RK3 (3rd order Implicit-Explicit)
    IMEXRK3,
}

/// Configuration for time stepping strategy
#[derive(Debug, Clone, Copy)]
pub struct SteppingConfig {
    /// Whether to use adaptive time stepping
    pub adaptive: bool,
    /// Whether to use local time stepping
    pub local: bool,
    /// Whether to use CFL-based time stepping
    pub use_cfl: bool,
}

impl Default for SteppingConfig {
    fn default() -> Self {
        Self {
            adaptive: true,
            local: false,
            use_cfl: true,
        }
    }
}

/// Parameters for time integration
#[derive(Debug, Clone)]
pub struct TimeIntegrationParams {
    /// Type of time integration method
    pub method: TimeIntegration,
    /// Time step size (if fixed)
    pub dt: Option<f64>,
    /// Final time
    pub t_final: f64,
    /// Maximum number of time steps
    pub max_steps: usize,
    /// Relative tolerance for adaptive time stepping
    pub rtol: f64,
    /// Absolute tolerance for adaptive time stepping
    pub atol: f64,
    /// Safety factor for adaptive time stepping
    pub safety_factor: f64,
    /// Minimum time step size
    pub dt_min: f64,
    /// Maximum time step size
    pub dt_max: f64,
    /// Stepping strategy configuration
    pub stepping: SteppingConfig,
    /// CFL number
    pub cfl: f64,
    /// Whether to output solution at each time step
    pub verbose: bool,
    /// Output interval (in time steps)
    pub output_interval: usize,
}

impl Default for TimeIntegrationParams {
    fn default() -> Self {
        Self {
            method: TimeIntegration::SSPRK3,
            dt: None,
            t_final: 1.0,
            max_steps: 1000,
            rtol: 1e-4,
            atol: 1e-6,
            safety_factor: 0.9,
            dt_min: 1e-10,
            dt_max: 0.1,
            stepping: SteppingConfig::default(),
            cfl: 0.1,
            verbose: false,
            output_interval: 10,
        }
    }
}

impl TimeIntegrationParams {
    /// Create a new set of time integration parameters
    pub fn new(method: TimeIntegration) -> Self {
        Self {
            method,
            ..Default::default()
        }
    }
    
    /// Set the time step size
    pub fn with_dt(mut self, dt: f64) -> Self {
        self.dt = Some(dt);
        self
    }
    
    /// Set the final time
    pub fn with_t_final(mut self, t_final: f64) -> Self {
        self.t_final = t_final;
        self
    }
    
    /// Set the maximum number of time steps
    pub fn with_max_steps(mut self, max_steps: usize) -> Self {
        self.max_steps = max_steps;
        self
    }
    
    /// Set the tolerance for adaptive time stepping
    pub fn with_tolerance(mut self, rtol: f64, atol: f64) -> Self {
        self.rtol = rtol;
        self.atol = atol;
        self
    }
    
    /// Set the safety factor for adaptive time stepping
    pub fn with_safety_factor(mut self, safety_factor: f64) -> Self {
        self.safety_factor = safety_factor;
        self
    }
    
    /// Set the minimum and maximum time step sizes
    pub fn with_dt_limits(mut self, dt_min: f64, dt_max: f64) -> Self {
        self.dt_min = dt_min;
        self.dt_max = dt_max;
        self
    }
    
    /// Enable or disable adaptive time stepping
    pub fn with_adaptive(mut self, adaptive: bool) -> Self {
        self.stepping.adaptive = adaptive;
        self
    }
    
    /// Enable or disable local time stepping
    pub fn with_local_time_stepping(mut self, local_time_stepping: bool) -> Self {
        self.stepping.local = local_time_stepping;
        self
    }
    
    /// Enable or disable CFL-based time stepping
    pub fn with_cfl(mut self, cfl: f64) -> Self {
        self.stepping.use_cfl = true;
        self.cfl = cfl;
        self
    }
    
    /// Enable or disable verbose output
    pub fn with_verbose(mut self, verbose: bool) -> Self {
        self.verbose = verbose;
        self
    }
    
    /// Set the output interval
    pub fn with_output_interval(mut self, output_interval: usize) -> Self {
        self.output_interval = output_interval;
        self
    }
}

/// Result of a time step
#[derive(Debug, Clone)]
pub struct TimeStepResult {
    /// Time at the end of the step
    pub t: f64,
    /// Time step size used
    pub dt: f64,
    /// Whether the step was accepted
    pub accepted: bool,
    /// Error estimate (for adaptive time stepping)
    pub error: Option<f64>,
    /// Number of function evaluations
    pub nfev: usize,
    /// Number of Jacobian evaluations
    pub njev: usize,
    /// Number of linear solves
    pub nlinsolve: usize,
    /// Whether the solver converged
    pub converged: bool,
    /// Error message (if any)
    pub message: Option<String>,
}

impl Default for TimeStepResult {
    fn default() -> Self {
        Self {
            t: 0.0,
            dt: 0.0,
            accepted: false,
            error: None,
            nfev: 0,
            njev: 0,
            nlinsolve: 0,
            converged: false,
            message: None,
        }
    }
}

/// Trait for time integration methods
pub trait TimeIntegrator: Send + Sync {
    /// Take a single time step
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)>;
    
    /// Get the order of the method
    fn order(&self) -> usize;
    
    /// Whether the method is implicit
    fn is_implicit(&self) -> bool;
    
    /// Whether the method is adaptive
    fn is_adaptive(&self) -> bool;
}

/// Forward Euler method (1st order explicit)
pub struct ForwardEuler;

impl TimeIntegrator for ForwardEuler {
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        _jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)>
    {
        let k1 = f(t, y)?;
        let y_new = y + k1 * dt;
        Ok((y_new, None))
    }
    
    fn order(&self) -> usize { 1 }
    fn is_implicit(&self) -> bool { false }
    fn is_adaptive(&self) -> bool { false }
}

/// Classic 4th order Runge-Kutta method
pub struct RK4;

impl TimeIntegrator for RK4 {
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        _jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)>
    {
        let k1 = f(t, y)?;
        let k2 = f(t + 0.5 * dt, &(y + &(&k1 * (0.5 * dt))))?;
        let k3 = f(t + 0.5 * dt, &(y + &(&k2 * (0.5 * dt))))?;
        let k4 = f(t + dt, &(y + &(&k3 * dt)))?;
        
        let y_new = y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        Ok((y_new, None))
    }
    
    fn order(&self) -> usize { 4 }
    fn is_implicit(&self) -> bool { false }
    fn is_adaptive(&self) -> bool { false }
}

/// Strong Stability Preserving RK3 (3rd order)
pub struct SSPRK3;

impl TimeIntegrator for SSPRK3 {
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        _jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)>
    {
        // Standard SSP-RK3 (Shu-Osher form)
        // u1 = u^n + dt * f(u^n)
        let k1 = f(t, y)?;
        let u1 = y + &(k1 * dt);
        
        // u2 = 3/4 * u^n + 1/4 * u1 + 1/4 * dt * f(u1)
        let k2 = f(t + dt, &u1)?;
        let u2 = (y * 0.75) + (u1 * 0.25) + (k2 * (0.25 * dt));
        
        // u_new = 1/3 * u^n + 2/3 * u2 + 2/3 * dt * f(u2)
        let k3 = f(t + 0.5 * dt, &u2)?;
        let y_new = (y * (1.0 / 3.0)) + (u2 * (2.0 / 3.0)) + (k3 * (2.0 / 3.0 * dt));
        
        Ok((y_new, None))
    }
    
    fn order(&self) -> usize { 3 }
    fn is_implicit(&self) -> bool { false }
    fn is_adaptive(&self) -> bool { false }
}

/// Implicit Euler method (1st order)
pub struct ImplicitEuler {
    /// Tolerance for Newton's method
    tol: f64,
    /// Maximum number of Newton iterations
    max_iter: usize,
}

impl Default for ImplicitEuler {
    fn default() -> Self {
        Self {
            tol: 1e-8,
            max_iter: 50,
        }
    }
}

impl ImplicitEuler {
    /// Create a new ImplicitEuler solver with custom parameters
    pub fn new(tol: f64, max_iter: usize) -> Self {
        Self { tol, max_iter }
    }
}

impl TimeIntegrator for ImplicitEuler {
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)>
    {
        // Initial guess (explicit Euler step)
        let f0 = f(t, y)?;
        let mut y_new = y + &(f0 * dt);
        
        // Newton iteration
        let mut iter = 0;
        let mut converged = false;
        
        while !converged && iter < self.max_iter {
            // Compute residual: G(y) = y - y_n - dt * f(t_{n+1}, y)
            let res = &y_new - y - &(f(t + dt, &y_new)? * dt);
            
            // Check convergence
            let res_norm = res.norm();
            if res_norm < self.tol {
                converged = true;
                break;
            }
            
            // Compute Jacobian: I - dt * J, where J = df/dy
            let jacobian = if let Some(jac_fn) = jac {
                jac_fn(t + dt, &y_new)?
            } else {
                // If no analytical Jacobian is provided, use finite differences
                self.finite_difference_jacobian(t + dt, &y_new, f)?
            };
            
            // Solve (I - dt * J) * dy = -res
            let jac_eye = DMatrix::identity(y.nrows(), y.ncols()) - &(jacobian * dt);
            let dy = jac_eye.lu().solve(&(-res.clone())).unwrap_or_else(|| -res);
            
            // Update solution
            y_new += &dy;
            iter += 1;
        }
        
        if !converged {
            return Err(DGError::NumericalError(
                format!("Newton iteration did not converge after {} iterations", self.max_iter)
            ));
        }
        
        Ok((y_new, None))
    }
    
    fn order(&self) -> usize { 1 }
    fn is_implicit(&self) -> bool { true }
    fn is_adaptive(&self) -> bool { false }
}

impl ImplicitEuler {
    /// Compute the Jacobian using finite differences
    fn finite_difference_jacobian(
        &self,
        t: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
    ) -> Result<DMatrix<f64>>
    {
        let eps = 1e-8;
        let n = y.len();
        let mut jac = DMatrix::zeros(n, n);
        let f0 = f(t, y)?;
        
        let mut y_pert = y.clone();
        for j in 0..n {
            // Perturb j-th component
            let yj = y[j];
            let h = eps * (1.0 + yj.abs());
            y_pert[j] = yj + h;
            
            // Compute perturbed function value
            let f_pert = f(t, &y_pert)?;
            
            // Finite difference approximation of the j-th column of the Jacobian
            for i in 0..n {
                jac[(i, j)] = (f_pert[i] - f0[i]) / h;
            }
            
            // Reset perturbation
            y_pert[j] = yj;
        }
        
        Ok(jac)
    }
}

/// Factory for creating time integrators
pub struct TimeIntegratorFactory;

impl TimeIntegratorFactory {
    /// Create a new time integrator
    pub fn create(method: TimeIntegration) -> Box<dyn TimeIntegrator> {
        match method {
            TimeIntegration::ForwardEuler => Box::new(ForwardEuler),
            TimeIntegration::RK4 => Box::new(RK4),
            TimeIntegration::ImplicitEuler => Box::new(ImplicitEuler::default()),
            _ => Box::new(SSPRK3), // Default to SSPRK3
        }
    }
}

/// Solver for time-dependent PDEs using DG methods
#[allow(dead_code)]
pub struct DGSolver {
    /// DG operator
    pub dg_op: DGOperator,
    /// Time integrator
    pub integrator: Box<dyn TimeIntegrator>,
    /// Time integration parameters
    pub params: TimeIntegrationParams,
    /// Current time
    pub t: f64,
    /// Current time step
    pub dt: f64,
    /// Current solution
    pub u: DMatrix<f64>,
    /// Number of time steps taken
    pub step_count: usize,
    /// Number of function evaluations
    pub nfev: usize,
    /// Number of Jacobian evaluations
    pub njev: usize,
    /// Number of linear solves
    pub nlinsolve: usize,
    /// Number of rejected steps
    pub nrejected: usize,
    /// Whether the solver is initialized
    pub initialized: bool,
}

#[allow(dead_code)]
impl DGSolver {
    /// Create a new DG solver
    pub fn new(
        dg_op: DGOperator,
        integrator: Box<dyn TimeIntegrator>,
        params: TimeIntegrationParams,
    ) -> Self {
        Self {
            dg_op,
            integrator,
            params,
            t: 0.0,
            dt: 0.0,
            u: DMatrix::zeros(0, 0),
            step_count: 0,
            nfev: 0,
            njev: 0,
            nlinsolve: 0,
            nrejected: 0,
            initialized: false,
        }
    }
    
    /// Initialize the solver with an initial condition
    pub fn initialize<F>(&mut self, u0: F) -> Result<()>
    where
        F: Fn(f64) -> DVector<f64>,
    {
        // Project the initial condition onto the DG basis
        self.u = self.dg_op.project(u0);
        
        // Set initial time step
        self.dt = self.compute_initial_dt()?;
        
        self.t = 0.0;
        self.step_count = 0;
        self.nfev = 0;
        self.njev = 0;
        self.nlinsolve = 0;
        self.nrejected = 0;
        self.initialized = true;
        
        Ok(())
    }
    
    /// Compute the initial time step
    fn compute_initial_dt(&self) -> Result<f64> {
        if let Some(dt) = self.params.dt {
            return Ok(dt);
        }
        
        if self.params.stepping.use_cfl {
            // Compute the maximum wave speed
            let max_wave_speed = self.compute_max_wave_speed()?;
            
            // Compute the CFL-based time step
            let h_min = 2.0 / (self.dg_op.order as f64 * self.dg_op.order as f64);
            let dt = self.params.cfl * h_min / (max_wave_speed + f64::EPSILON);
            
            Ok(dt.min(self.params.dt_max).max(self.params.dt_min))
        } else {
            // Default to a small time step
            Ok(self.params.dt_min)
        }
    }
    
    /// Compute the maximum wave speed in the domain
    fn compute_max_wave_speed(&self) -> Result<f64> {
        // This is a simplified version that assumes a constant wave speed
        // In practice, this should be computed based on the solution
        Ok(1.0)
    }
    
    /// Take a single time step
    pub fn step<F, J>(
        &mut self,
        f: &F,
        jac: Option<&J>,
    ) -> Result<TimeStepResult>
    where
        F: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
        J: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
    {
        if !self.initialized {
            return Err(DGError::NumericalError(
                "Solver not initialized. Call initialize() first.".to_string()
            ));
        }
        
        if self.t >= self.params.t_final {
            return Ok(TimeStepResult {
                t: self.t,
                dt: 0.0,
                accepted: false,
                error: None,
                nfev: 0,
                njev: 0,
                nlinsolve: 0,
                converged: true,
                message: Some("Simulation completed successfully".to_string()),
            });
        }
        
        let _t_start = Instant::now();
        
        // Adjust time step to hit t_final exactly
        let dt = if self.t + self.dt > self.params.t_final {
            self.params.t_final - self.t
        } else {
            self.dt
        };
        
        // Take a time step
        let jac_dyn = jac.map(|j| j as &dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>);
        let (u_new, error) = self.integrator.step(self.t, dt, &self.u, f, jac_dyn)?;
        
        // Update statistics
        self.nfev += 1; // This is a simplification; actual count depends on the method
        
        // Accept the step
        self.u = u_new;
        self.t += dt;
        self.step_count += 1;
        
        // Update time step for next iteration
        if self.params.stepping.adaptive {
            if let Some(err) = error {
                // Compute the optimal time step
                let scale = (self.params.safety_factor * self.params.rtol / (err + f64::EPSILON))
                    .powf(1.0 / (self.integrator.order() as f64 + 1.0));
                
                self.dt = (dt * scale)
                    .max(self.params.dt_min)
                    .min(self.params.dt_max);
            }
        }
        
        // Output progress
        if self.params.verbose && (self.step_count.is_multiple_of(self.params.output_interval) || self.t >= self.params.t_final) {
            println!(
                "Step {}: t = {:.6}, dt = {:.2e}, |u| = {:.6e}",
                self.step_count,
                self.t,
                dt,
                self.u.norm()
            );
        }
        
        Ok(TimeStepResult {
            t: self.t,
            dt,
            accepted: true,
            error,
            nfev: 1, // Simplified
            njev: 0,  // Simplified
            nlinsolve: 0, // Simplified
            converged: true,
            message: None,
        })
    }
    
    /// Run the solver until t_final is reached
    pub fn solve<F, J>(
        &mut self,
        f: F,
        jac: Option<J>,
    ) -> Result<()>
    where
        F: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
        J: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
    {
        if !self.initialized {
            return Err(DGError::NumericalError(
                "Solver not initialized. Call initialize() first.".to_string()
            ));
        }
        
        let start_time = Instant::now();
        
        if self.params.verbose {
            println!("Starting time integration...");
            println!("  Method: {:?}", self.params.method);
            println!("  t_final: {}", self.params.t_final);
            println!("  Initial dt: {}", self.dt);
            println!("  Max steps: {}", self.params.max_steps);
        }
        
        // Main time-stepping loop
        let mut result = TimeStepResult::default();
        
        while self.t < self.params.t_final && self.step_count < self.params.max_steps {
            result = self.step(&f, jac.as_ref())?;
            
            if !result.converged {
                if self.params.verbose {
                    println!("Warning: Time step failed at t = {}", self.t);
                }
                
                if self.dt <= self.params.dt_min {
                    return Err(DGError::NumericalError(
                        format!("Time step too small at t = {}", self.t)
                    ));
                }
                
                // Reduce time step and try again
                self.dt = (self.dt / 2.0).max(self.params.dt_min);
                self.nrejected += 1;
                
                if self.params.verbose {
                    println!("  Reducing time step to dt = {}", self.dt);
                }
                
                continue;
            }
            
            self.step_count += 1;
        }
        
        let elapsed = start_time.elapsed();
        
        if self.params.verbose {
            println!("Time integration completed in {:.3} seconds", elapsed.as_secs_f64());
            println!("  Final time: {}", self.t);
            println!("  Time steps: {}", self.step_count);
            println!("  Rejected steps: {}", self.nrejected);
            println!("  Function evaluations: {}", self.nfev);
            println!("  Jacobian evaluations: {}", self.njev);
            println!("  Linear solves: {}", self.nlinsolve);
        }
        
        if self.t < self.params.t_final && self.step_count >= self.params.max_steps {
            return Err(DGError::NumericalError(
                format!("Maximum number of time steps ({}) reached", self.params.max_steps)
            ));
        }
        
        Ok(())
    }
    
    /// Get the current solution
    pub fn solution(&self) -> &DMatrix<f64> {
        &self.u
    }
    
    /// Evaluate the solution at a point x in the reference element [-1, 1]
    pub fn evaluate(&self, x: f64) -> DVector<f64> {
        let num_components = self.u.nrows();
        let num_basis = self.u.ncols();
        let mut u = DVector::zeros(num_components);
        
        for i in 0..num_basis {
            let phi_i = legendre_poly(i, x);
            for c in 0..num_components {
                u[c] += self.u[(c, i)] * phi_i;
            }
        }
        
        u
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    use crate::high_order::{DGOperatorParams, FluxType, LimiterType};
    
    #[test]
    fn test_forward_euler() {
        // Test the Forward Euler method on the ODE: du/dt = -u, u(0) = 1
        // Exact solution: u(t) = exp(-t)
        
        let f = |_t: f64, u: &DMatrix<f64>| {
            Ok(-u.clone())
        };
        
        let dt = 0.01;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;
        
        let mut u = DMatrix::from_element(1, 1, 1.0); // u(0) = 1
        
        let integrator = ForwardEuler;
        
        let jac: Option<&dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-2);
    }

    #[test]
    fn test_rk4() {
        // Test the RK4 method on the ODE: du/dt = -u, u(0) = 1
        // Exact solution: u(t) = exp(-t)

        let f = |_t: f64, u: &DMatrix<f64>| {
            Ok(-u.clone())
        };

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0); // u(0) = 1

        let integrator = RK4;

        let jac: Option<&dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-6);
    }

    #[test]
    fn test_ssprk3() {
        // Test the SSPRK3 method on the ODE: du/dt = -u, u(0) = 1
        // Exact solution: u(t) = exp(-t)

        let f = |_t: f64, u: &DMatrix<f64>| {
            Ok(-u.clone())
        };

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0); // u(0) = 1

        let integrator = SSPRK3;

        let jac: Option<&dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }
        
        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-4);
    }
    
    #[test]
    fn test_implicit_euler() {
        // Test the Implicit Euler method on the ODE: du/dt = -u, u(0) = 1
        // Exact solution: u(t) = exp(-t)
        
        let f = |_t: f64, u: &DMatrix<f64>| {
            Ok(-u.clone())
        };
        
        // For the linear problem du/dt = A*u, the Jacobian is A
        let jac = |_t: f64, _u: &DMatrix<f64>| {
            Ok(DMatrix::from_element(1, 1, -1.0))
        };
        
        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;
        
        let mut u = DMatrix::from_element(1, 1, 1.0); // u(0) = 1
        
        let integrator = ImplicitEuler::default();
        
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, Some(&jac)).unwrap();
            u = u_new;
        }
        
        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 2e-2);
    }
    
    #[test]
    fn test_dg_solver() -> Result<()> {
        // Test the DG solver on the ODE: du/dt = -u
        // with initial condition u(x,0) = 1 + x + x^2
        // Exact solution: u(x,t) = (1 + x + x^2) * exp(-t)
        
        // Create a DG operator
        let order = 2;
        let num_components = 1;
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::Central)
            .with_surface_flux(FluxType::LaxFriedrichs)
            .with_limiter(LimiterType::None);
        
        let dg_op = DGOperator::new(order, num_components, Some(params))?;
        
        // Create a time integrator
        let integrator = TimeIntegratorFactory::create(TimeIntegration::SSPRK3);
        
        // Set up the solver
        let t_final = 1.0;
        let solver_params = TimeIntegrationParams::new(TimeIntegration::SSPRK3)
            .with_t_final(t_final)
            .with_dt(0.01)
            .with_verbose(false);
        
        let mut solver = DGSolver::new(dg_op, integrator, solver_params);
        
        // Initial condition: u(x,0) = 1 + x + x^2
        let u0 = |x: f64| DVector::from_vec(vec![1.0 + x + x * x]);
        
        // Initialize the solver
        solver.initialize(u0)?;
        
        // Define the right-hand side function: du/dt = -u
        let f = |_t: f64, u: &DMatrix<f64>| {
            Ok(-u.clone())
        };
        
        // Run the solver (no Jacobian needed for explicit methods)
        solver.solve(f, None::<fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>>)?;
        
        // Check the solution at the final time
        let x = 0.5; // Test point
        let u_num = solver.evaluate(x)[0];
        let u_exact = (1.0 + x + x * x) * (-t_final).exp();
        
        assert_relative_eq!(u_num, u_exact, epsilon = 1e-3);
        
        Ok(())
    }
}
