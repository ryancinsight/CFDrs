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

/// Minimum time step threshold to avoid numerical issues
const MIN_DT_THRESHOLD: f64 = 1e-10;

/// Factor for percentage conversion
const PERCENTAGE_FACTOR: f64 = 100.0;

/// State for PISO solver execution
pub struct PisoState<T: RealField + Copy> {
    /// Convergence monitor
    pub monitor: ConvergenceMonitor<T>,
    /// Buffer for double-buffering pattern (allocated once)
    pub fields_buffer: SimulationFields<T>,
}

impl<T: RealField + Copy> PisoState<T> {
    /// Create new state with initialized fields
    #[must_use]
    pub fn new(fields: &SimulationFields<T>) -> Self {
        Self {
            monitor: ConvergenceMonitor::new(),
            fields_buffer: fields.clone(),
        }
    }
}

/// PISO solver for incompressible flow (transient algorithm)
/// Stateless solver - all mutable state is external
pub struct PisoSolver<T: RealField + Copy> {
    /// Solver configuration
    config: PisoConfig<T>,
    /// Velocity predictor
    predictor: VelocityPredictor<T>,
    /// Pressure corrector
    corrector: PressureCorrector<T>,
    /// Convergence criteria (for inner iterations within a time step)
    criteria: ConvergenceCriteria<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::iter::Sum> PisoSolver<T> {
    /// Create new PISO solver
    pub fn new(config: PisoConfig<T>, grid: &StructuredGrid2D<T>) -> Self {
        let predictor = VelocityPredictor::new(grid, config.velocity_relaxation);
        let corrector =
            PressureCorrector::new(grid, config.n_correctors, config.pressure_relaxation);
        let criteria = ConvergenceCriteria::default();

        Self {
            config,
            predictor,
            corrector,
            criteria,
        }
    }

    /// Advance solution by one time step
    /// This is the core PISO algorithm: predictor + corrector(s) for a single time step
    pub fn advance_one_step(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
    ) -> Result<()> {
        self.advance_with_dt(fields, grid, state, self.config.time_step)
    }

    /// Advance solution by a specified time step
    pub fn advance_with_dt(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
        dt: T,
    ) -> Result<()> {
        // Copy current state to buffer for residual calculation
        // This is much cheaper than cloning in every iteration
        state.fields_buffer.copy_from(fields)?;

        // Step 1: Velocity predictor
        self.predictor.predict(fields, dt)?;

        // Step 2: Pressure correction (multiple correctors)
        // The corrector internally handles the number of correction steps
        self.corrector.correct(fields, dt)?;

        // Step 3: Update convergence monitor for diagnostics
        // Note: This is for monitoring within-timestep convergence, not steady-state
        state
            .monitor
            .update(&state.fields_buffer, fields, grid.nx, grid.ny);
        state.monitor.iteration += 1;

        Ok(())
    }

    /// Run transient simulation for specified number of time steps with callback
    pub fn solve_transient_with_callback<F>(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
        num_steps: usize,
        mut on_step_complete: F,
    ) -> Result<()>
    where
        F: FnMut(usize, &SimulationFields<T>) -> Result<()>,
    {
        for step in 0..num_steps {
            self.advance_one_step(fields, grid, state)?;

            // Optional: Log progress based on configuration
            if let Some(freq) = self.config.log_frequency {
                if freq > 0 && step % freq == 0 && step > 0 {
                    if let Some(vel_res) = state.monitor.velocity_residuals.last() {
                        // Convert to f64 for display since T might not implement LowerExp
                        if let Some(vel_res_f64) = vel_res.to_f64() {
                            tracing::info!("Step {}: velocity residual = {:e}", step, vel_res_f64);
                        }
                    }
                }
            }

            // Call user-provided callback
            on_step_complete(step, fields)?;
        }
        Ok(())
    }

    /// Run transient simulation for specified number of time steps (backward compatible)
    pub fn solve_transient(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
        num_steps: usize,
    ) -> Result<()> {
        self.solve_transient_with_callback(fields, grid, state, num_steps, |_, _| Ok(()))
    }

    /// Run transient simulation for specified physical duration
    pub fn solve_for_duration(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
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
            let min_dt_threshold = T::from_f64(MIN_DT_THRESHOLD)
                .expect("Failed to represent minimum dt threshold in numeric type T");
            if dt <= min_dt_threshold {
                break;
            }

            // Advance by the calculated time step
            self.advance_with_dt(fields, grid, state, dt)?;
            current_time += dt;
            step += 1;

            // Optional: Log progress based on configuration
            if let Some(freq) = self.config.log_frequency {
                if freq > 0 && step % freq == 0 {
                    let percent_factor = T::from_f64(PERCENTAGE_FACTOR)
                        .expect("Failed to represent percentage factor in numeric type T");
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
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        state: &mut PisoState<T>,
        max_steps: usize,
    ) -> Result<bool> {
        for step in 0..max_steps {
            self.advance_one_step(fields, grid, state)?;

            // Check if solution has reached steady state
            // (residuals below tolerance)
            if state.monitor.is_converged(&self.criteria) {
                tracing::info!("Reached steady state at step {}", step);
                return Ok(true);
            }
        }

        tracing::warn!("Max steps {} reached without steady state", max_steps);
        Ok(false)
    }

    /// Reset convergence history
    pub fn reset_history(&mut self, state: &mut PisoState<T>) {
        state.monitor = ConvergenceMonitor::new();
    }

    /// Get convergence history
    pub fn convergence_history<'a>(&self, state: &'a PisoState<T>) -> &'a ConvergenceMonitor<T> {
        &state.monitor
    }

    /// Update configuration and clear buffer if needed
    pub fn set_config(&mut self, config: PisoConfig<T>) {
        self.config = config;
    }

    /// Get current configuration
    pub fn config(&self) -> &PisoConfig<T> {
        &self.config
    }
}
