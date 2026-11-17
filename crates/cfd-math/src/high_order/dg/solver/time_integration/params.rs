//! Configuration and result types for DG time integration.
//!
//! Split out of `time_integration.rs` so the integrator implementations keep
//! their own file; every public path is re-exported from the parent module.

/// Type of time integration method
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TimeIntegration {
    /// Forward Euler (1st order explicit)
    ForwardEuler,
    /// Classic 4th order Runge-Kutta
    RK4,
    /// Strong Stability Preserving RK3 (3rd order)
    SSPRK3,
    /// Implicit Euler (1st order)
    ImplicitEuler,
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
