//! Main PISO solver implementation

use super::{
    config::PisoConfig,
    convergence::{ConvergenceCriteria, ConvergenceMonitor},
    corrector::PressureCorrector,
    predictor::VelocityPredictor,
};
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

/// PISO solver for incompressible flow (transient algorithm)
pub struct PisoSolver<T: RealField + Copy> {
    /// Solver configuration
    config: PisoConfig<T>,
    /// Velocity predictor
    predictor: VelocityPredictor<T>,
    /// Pressure corrector
    corrector: PressureCorrector<T>,
    /// Convergence monitor
    monitor: ConvergenceMonitor<T>,
    /// Convergence criteria (for inner iterations within a time step)
    criteria: ConvergenceCriteria<T>,
    /// Buffer for double-buffering pattern (allocated once)
    fields_buffer: Option<SimulationFields<T>>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::iter::Sum> PisoSolver<T> {
    /// Create new PISO solver
    pub fn new(config: PisoConfig<T>, grid: &StructuredGrid2D<T>) -> Self {
        let predictor = VelocityPredictor::new(grid, config.velocity_relaxation);
        let corrector =
            PressureCorrector::new(grid, config.n_correctors, config.pressure_relaxation);
        let monitor = ConvergenceMonitor::new();
        let criteria = ConvergenceCriteria::default();

        Self {
            config,
            predictor,
            corrector,
            monitor,
            criteria,
            fields_buffer: None,
        }
    }

    /// Advance solution by one time step
    /// This is the core PISO algorithm: predictor + corrector(s) for a single time step
    pub fn advance_one_step(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<()> {
        self.advance_with_dt(fields, grid, self.config.time_step)
    }

    /// Advance solution by a specified time step
    pub fn advance_with_dt(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        dt: T,
    ) -> Result<()> {
        // Initialize buffer on first use (lazy allocation)
        if self.fields_buffer.is_none() {
            self.fields_buffer = Some(fields.clone());
        }

        // Get reference to buffer
        let fields_old = self.fields_buffer.as_mut().unwrap();

        // Copy current state to buffer for residual calculation
        // This is much cheaper than cloning in every iteration
        fields_old.copy_from(fields);

        // Step 1: Velocity predictor
        self.predictor.predict(fields, dt)?;

        // Step 2: Pressure correction (multiple correctors)
        // The corrector internally handles the number of correction steps
        self.corrector.correct(fields, dt)?;

        // Step 3: Update convergence monitor for diagnostics
        // Note: This is for monitoring within-timestep convergence, not steady-state
        self.monitor.update(fields_old, fields, grid.nx, grid.ny);
        self.monitor.iteration += 1;

        Ok(())
    }

    /// Run transient simulation for specified number of time steps with callback
    pub fn solve_transient_with_callback<F>(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        num_steps: usize,
        mut on_step_complete: F,
    ) -> Result<()>
    where
        F: FnMut(usize, &SimulationFields<T>) -> Result<()>,
    {
        for step in 0..num_steps {
            self.advance_one_step(fields, grid)?;

            // Optional: Log progress based on configuration
            if let Some(freq) = self.config.log_frequency {
                if freq > 0 && step % freq == 0 && step > 0 {
                    if let Some(vel_res) = self.monitor.velocity_residuals.last() {
                        // Convert to f64 for display since T might not implement LowerExp
                        if let Some(vel_res_f64) = vel_res.to_f64() {
                            tracing::info!("Step {}: velocity residual = {:e}", step, vel_res_f64);
                        }
                    }
                }
            }

            // Call user-provided callback
            on_step_complete(step, fields)?
        }
        Ok(())
    }

    /// Run transient simulation for specified number of time steps (backward compatible)
    pub fn solve_transient(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        num_steps: usize,
    ) -> Result<()> {
        self.solve_transient_with_callback(fields, grid, num_steps, |_, _| Ok(()))
    }

    /// Run transient simulation for specified physical duration
    pub fn solve_for_duration(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        total_duration: T,
    ) -> Result<()> {
        let mut current_time = T::zero();
        let mut step = 0;

        while current_time < total_duration {
            // Check if adding another time_step would overshoot, and adjust if needed
            let remaining = total_duration - current_time;
            let dt = if remaining < self.config.time_step {
                remaining
            } else {
                self.config.time_step
            };

            // Skip if time step becomes too small
            let min_dt_threshold = T::from_f64(1e-10)
                .expect("Failed to represent minimum dt threshold (1e-10) in numeric type T");
            if dt <= min_dt_threshold {
                break;
            }

            // Advance by the calculated time step
            self.advance_with_dt(fields, grid, dt)?;
            current_time = current_time + dt;
            step += 1;

            // Optional: Log progress based on configuration
            if let Some(freq) = self.config.log_frequency {
                if freq > 0 && step % freq == 0 {
                    let percent_factor = T::from_f64(100.0)
                        .expect("Failed to represent percentage factor (100.0) in numeric type T");
                    let progress = current_time / total_duration * percent_factor;
                    // Convert to f64 for display
                    if let (Some(progress_f64), Some(time_f64)) =
                        (progress.to_f64(), current_time.to_f64())
                    {
                        tracing::info!(
                            "Simulation progress: {:.1}% (t = {:.3})",
                            progress_f64,
                            time_f64
                        );
                    }
                }
            }
        }

        Ok(())
    }

    /// Run simulation until steady state is reached (if applicable)
    /// Note: PISO is primarily for transient problems, but this can be useful
    /// for finding steady-state solutions of certain problems
    pub fn solve_to_steady_state(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        max_steps: usize,
    ) -> Result<bool> {
        for step in 0..max_steps {
            self.advance_one_step(fields, grid)?;

            // Check if solution has reached steady state
            // (residuals below tolerance)
            if self.monitor.is_converged(&self.criteria) {
                tracing::info!("Reached steady state at step {}", step);
                return Ok(true);
            }
        }

        tracing::warn!("Max steps {} reached without steady state", max_steps);
        Ok(false)
    }

    /// Reset convergence history
    pub fn reset_history(&mut self) {
        self.monitor = ConvergenceMonitor::new();
    }

    /// Get convergence history
    pub fn convergence_history(&self) -> &ConvergenceMonitor<T> {
        &self.monitor
    }

    /// Update configuration
    pub fn set_config(&mut self, config: PisoConfig<T>) {
        self.config = config;
        // Clear buffer if grid size might have changed
        self.fields_buffer = None;
    }

    /// Get current configuration
    pub fn config(&self) -> &PisoConfig<T> {
        &self.config
    }
}
