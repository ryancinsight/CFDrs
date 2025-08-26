//! Main PISO solver implementation

use super::{
    config::PisoConfig,
    convergence::{ConvergenceCriteria, ConvergenceMonitor},
    corrector::PressureCorrector,
    predictor::VelocityPredictor,
};
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use cfd_core::Result;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// PISO solver for incompressible flow
pub struct PisoSolver<T: RealField + Copy> {
    /// Solver configuration
    config: PisoConfig<T>,
    /// Velocity predictor
    predictor: VelocityPredictor<T>,
    /// Pressure corrector
    corrector: PressureCorrector<T>,
    /// Convergence monitor
    monitor: ConvergenceMonitor<T>,
    /// Convergence criteria
    criteria: ConvergenceCriteria<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy + std::iter::Sum> PisoSolver<T> {
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
        }
    }

    /// Perform one PISO iteration
    pub fn iterate(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<bool> {
        // Store old fields for convergence check
        let fields_old = fields.clone();

        // Step 1: Velocity predictor
        self.predictor.predict(fields, self.config.time_step)?;

        // Step 2: Pressure correction (multiple correctors)
        self.corrector.correct(fields, self.config.time_step)?;

        // Step 3: Update convergence monitor
        self.monitor.update(&fields_old, fields, grid.nx, grid.ny);

        // Check convergence
        Ok(self.monitor.is_converged(&self.criteria))
    }

    /// Run PISO solver until convergence
    pub fn solve(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<()> {
        self.monitor = ConvergenceMonitor::new();

        loop {
            let converged = self.iterate(fields, grid)?;

            if converged {
                break;
            }

            if self.monitor.iteration >= self.criteria.max_iterations {
                // Warning: reached max iterations without convergence
                break;
            }
        }

        Ok(())
    }

    /// Get convergence history
    pub fn convergence_history(&self) -> &ConvergenceMonitor<T> {
        &self.monitor
    }

    /// Update configuration
    pub fn set_config(&mut self, config: PisoConfig<T>) {
        self.config = config;
    }
}
