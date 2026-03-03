//! DG time-dependent PDE solver.
//!
//! Combines a [`DGOperator`] with a [`TimeIntegrator`] to evolve DG
//! solutions forward in time using configurable explicit or implicit schemes.

mod time_integration;

pub use time_integration::*;

use super::{legendre_poly, DGError, DGOperator, Result};
use nalgebra::{DMatrix, DVector};
use std::time::Instant;

/// Solver for time-dependent PDEs using DG methods
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

    /// Compute the maximum wave speed in the domain.
    ///
    /// Scans all element-interface pairs in the current solution to obtain
    /// a global upper bound on the characteristic speed, which controls the
    /// CFL-based time step.
    fn compute_max_wave_speed(&self) -> Result<f64> {
        use super::flux::{FluxFactory, NumericalFlux};

        let n_elem = self.u.ncols();
        let n_dof = self.u.nrows();
        if n_elem == 0 || n_dof == 0 {
            return Ok(f64::EPSILON);
        }

        let flux = FluxFactory::create(self.dg_op.params.surface_flux);
        // Unit normal (1-D reference direction)
        let n = DVector::from_element(1, 1.0);

        let mut max_speed: f64 = f64::EPSILON;

        for e in 0..n_elem.saturating_sub(1) {
            // Right boundary of element e
            let u_l = DVector::from_element(1, self.u[(n_dof - 1, e)]);
            // Left boundary of element e+1
            let u_r = DVector::from_element(1, self.u[(0, e + 1)]);
            let speed = flux.max_wave_speed(&u_l, &u_r, &n);
            if speed > max_speed {
                max_speed = speed;
            }
        }

        Ok(max_speed)
    }

    /// Take a single time step
    pub fn step<F, J>(&mut self, f: &F, jac: Option<&J>) -> Result<TimeStepResult>
    where
        F: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
        J: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
    {
        if !self.initialized {
            return Err(DGError::NumericalError(
                "Solver not initialized. Call initialize() first.".to_string(),
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
        let step_fev = self.integrator.stages();
        self.nfev += step_fev;

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

                self.dt = (dt * scale).max(self.params.dt_min).min(self.params.dt_max);
            }
        }

        // Output progress
        if self.params.verbose
            && (self.step_count.is_multiple_of(self.params.output_interval)
                || self.t >= self.params.t_final)
        {
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
            nfev: step_fev,
            njev: usize::from(self.integrator.is_implicit()),
            nlinsolve: usize::from(self.integrator.is_implicit()),
            converged: true,
            message: None,
        })
    }

    /// Run the solver until t_final is reached
    pub fn solve<F, J>(&mut self, f: F, jac: Option<J>) -> Result<()>
    where
        F: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
        J: Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>,
    {
        if !self.initialized {
            return Err(DGError::NumericalError(
                "Solver not initialized. Call initialize() first.".to_string(),
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
        while self.t < self.params.t_final && self.step_count < self.params.max_steps {
            let result = self.step(&f, jac.as_ref())?;

            if !result.converged {
                if self.params.verbose {
                    println!("Warning: Time step failed at t = {}", self.t);
                }

                if self.dt <= self.params.dt_min {
                    return Err(DGError::NumericalError(format!(
                        "Time step too small at t = {}",
                        self.t
                    )));
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
            println!(
                "Time integration completed in {:.3} seconds",
                elapsed.as_secs_f64()
            );
            println!("  Final time: {}", self.t);
            println!("  Time steps: {}", self.step_count);
            println!("  Rejected steps: {}", self.nrejected);
            println!("  Function evaluations: {}", self.nfev);
            println!("  Jacobian evaluations: {}", self.njev);
            println!("  Linear solves: {}", self.nlinsolve);
        }

        if self.t < self.params.t_final && self.step_count >= self.params.max_steps {
            return Err(DGError::NumericalError(format!(
                "Maximum number of time steps ({}) reached",
                self.params.max_steps
            )));
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
    type Jacobian = dyn Fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>;

    #[test]
    fn test_forward_euler() {
        let f = |_t: f64, u: &DMatrix<f64>| Ok(-u.clone());

        let dt = 0.01;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0);

        let integrator = ForwardEuler;

        let jac: Option<&Jacobian> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-2);
    }

    #[test]
    fn test_rk4() {
        let f = |_t: f64, u: &DMatrix<f64>| Ok(-u.clone());

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0);

        let integrator = RK4;

        let jac: Option<&Jacobian> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-6);
    }

    #[test]
    fn test_ssprk3() {
        let f = |_t: f64, u: &DMatrix<f64>| Ok(-u.clone());

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0);

        let integrator = SSPRK3;

        let jac: Option<&Jacobian> = None;
        for _ in 0..n_steps {
            let (u_new, _) = integrator.step(0.0, dt, &u, &f, jac).unwrap();
            u = u_new;
        }

        let exact = (-t_final).exp();
        assert_relative_eq!(u[(0, 0)], exact, epsilon = 1e-4);
    }

    #[test]
    fn test_implicit_euler() {
        let f = |_t: f64, u: &DMatrix<f64>| Ok(-u.clone());
        let jac = |_t: f64, _u: &DMatrix<f64>| Ok(DMatrix::from_element(1, 1, -1.0));

        let dt = 0.1;
        let t_final = 1.0;
        let n_steps = (t_final / dt) as usize;

        let mut u = DMatrix::from_element(1, 1, 1.0);

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
        let order = 2;
        let num_components = 1;
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::Central)
            .with_surface_flux(FluxType::LaxFriedrichs)
            .with_limiter(LimiterType::None);

        let dg_op = DGOperator::new(order, num_components, Some(params))?;

        let integrator = TimeIntegratorFactory::create(TimeIntegration::SSPRK3);

        let t_final = 1.0;
        let solver_params = TimeIntegrationParams::new(TimeIntegration::SSPRK3)
            .with_t_final(t_final)
            .with_dt(0.01)
            .with_verbose(false);

        let mut solver = DGSolver::new(dg_op, integrator, solver_params);

        let u0 = |x: f64| DVector::from_vec(vec![1.0 + x + x * x]);
        solver.initialize(u0)?;

        let f = |_t: f64, u: &DMatrix<f64>| Ok(-u.clone());
        solver.solve(f, None::<fn(f64, &DMatrix<f64>) -> Result<DMatrix<f64>>>)?;

        let x = 0.5;
        let u_num = solver.evaluate(x)[0];
        let u_exact = (1.0 + x + x * x) * (-t_final).exp();

        assert_relative_eq!(u_num, u_exact, epsilon = 1e-3);

        Ok(())
    }
}
