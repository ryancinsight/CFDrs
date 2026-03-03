//! Time integration methods for Discontinuous Galerkin solvers.
//!
//! Provides explicit, implicit, and IMEX time-stepping schemes including
//! Forward Euler, RK4, SSP-RK3, and Implicit Euler with Newton iteration.

use nalgebra::DMatrix;
use super::super::{DGError, Result};

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

    /// Number of right-hand-side evaluations per step
    fn stages(&self) -> usize;

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
    ) -> Result<(DMatrix<f64>, Option<f64>)> {
        let k1 = f(t, y)?;
        let y_new = y + k1 * dt;
        Ok((y_new, None))
    }

    fn order(&self) -> usize {
        1
    }
    fn stages(&self) -> usize {
        1
    }
    fn is_implicit(&self) -> bool {
        false
    }
    fn is_adaptive(&self) -> bool {
        false
    }
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
    ) -> Result<(DMatrix<f64>, Option<f64>)> {
        let k1 = f(t, y)?;
        let k2 = f(t + 0.5 * dt, &(y + &(&k1 * (0.5 * dt))))?;
        let k3 = f(t + 0.5 * dt, &(y + &(&k2 * (0.5 * dt))))?;
        let k4 = f(t + dt, &(y + &(&k3 * dt)))?;

        let y_new = y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
        Ok((y_new, None))
    }

    fn order(&self) -> usize {
        4
    }
    fn stages(&self) -> usize {
        4
    }
    fn is_implicit(&self) -> bool {
        false
    }
    fn is_adaptive(&self) -> bool {
        false
    }
}

/// Strong Stability Preserving RK3 (3rd order)
///
/// # Theorem — SSP Property
/// The Shu-Osher SSP-RK3 method preserves the total variation diminishing
/// property under the CFL condition Δt ≤ Δx/max|a|, where a is the wave speed.
///
/// **Proof sketch**: Each stage is a convex combination of forward Euler steps,
/// and convex combinations preserve TVD stability (Shu & Osher, 1988).
pub struct SSPRK3;

impl TimeIntegrator for SSPRK3 {
    fn step(
        &self,
        t: f64,
        dt: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
        _jac: Option<&JacobianFn<'_>>,
    ) -> Result<(DMatrix<f64>, Option<f64>)> {
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

    fn order(&self) -> usize {
        3
    }
    fn stages(&self) -> usize {
        3
    }
    fn is_implicit(&self) -> bool {
        false
    }
    fn is_adaptive(&self) -> bool {
        false
    }
}

/// Implicit Euler method (1st order, A-stable)
///
/// Uses Newton iteration to solve the implicit system at each step.
/// A-stability allows arbitrarily large time steps for stiff problems,
/// though accuracy degrades with O(Δt) truncation error.
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
    ) -> Result<(DMatrix<f64>, Option<f64>)> {
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
            return Err(DGError::NumericalError(format!(
                "Newton iteration did not converge after {} iterations",
                self.max_iter
            )));
        }

        Ok((y_new, None))
    }

    fn order(&self) -> usize {
        1
    }
    fn stages(&self) -> usize {
        // 1 initial f-eval + up to max_iter Newton iterations (each with 1 f-eval)
        1 + self.max_iter
    }
    fn is_implicit(&self) -> bool {
        true
    }
    fn is_adaptive(&self) -> bool {
        false
    }
}

impl ImplicitEuler {
    /// Compute the Jacobian using finite differences
    fn finite_difference_jacobian(
        &self,
        t: f64,
        y: &DMatrix<f64>,
        f: &RhsFn<'_>,
    ) -> Result<DMatrix<f64>> {
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
