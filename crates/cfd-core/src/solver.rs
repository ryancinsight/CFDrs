//! Solver trait and common solver functionality.

use crate::{problem::Problem, Result};
use nalgebra::RealField;

/// Main solver trait for CFD simulations
pub trait Solver<T: RealField>: Send + Sync {
    /// Problem type this solver can handle
    type Problem: Problem<T>;
    /// Solution type produced by this solver
    type Solution;

    /// Solve the given problem
    fn solve(&mut self, problem: &Self::Problem) -> Result<Self::Solution>;

    /// Get solver name
    fn name(&self) -> &str;

    /// Get solver configuration
    fn config(&self) -> &SolverConfig<T>;

    /// Set solver configuration
    fn set_config(&mut self, config: SolverConfig<T>);
}

/// Common solver configuration
#[derive(Debug, Clone)]
pub struct SolverConfig<T: RealField> {
    /// Convergence tolerance
    pub tolerance: T,
    /// Maximum number of iterations
    pub max_iterations: usize,
    /// Relaxation factor
    pub relaxation_factor: T,
    /// Verbosity level (0 = silent, 1 = summary, 2 = detailed)
    pub verbosity: u8,
    /// Use parallel execution
    pub parallel: bool,
    /// Number of threads (None = use all available)
    pub num_threads: Option<usize>,
}

impl<T: RealField> Default for SolverConfig<T> {
    fn default() -> Self {
        Self {
            tolerance: T::from(1e-6).unwrap(),
            max_iterations: 1000,
            relaxation_factor: T::one(),
            verbosity: 1,
            parallel: true,
            num_threads: None,
        }
    }
}

/// Iterative solver trait
pub trait IterativeSolver<T: RealField>: Solver<T> {
    /// Perform one iteration
    fn iterate(&mut self, state: &mut Self::Solution) -> Result<T>;

    /// Check convergence
    fn is_converged(&self, residual: T) -> bool {
        residual < self.config().tolerance
    }

    /// Run the iterative solver
    fn solve_iterative(&mut self, problem: &Self::Problem) -> Result<Self::Solution> {
        let mut state = self.initialize(problem)?;
        let config = self.config();
        
        for iter in 0..config.max_iterations {
            let residual = self.iterate(&mut state)?;
            
            if config.verbosity >= 2 {
                tracing::debug!("Iteration {}: residual = {:e}", iter + 1, residual);
            }
            
            if self.is_converged(residual) {
                if config.verbosity >= 1 {
                    tracing::info!(
                        "Converged after {} iterations (residual: {:e})",
                        iter + 1,
                        residual
                    );
                }
                return Ok(state);
            }
        }
        
        Err(crate::Error::ConvergenceFailure {
            iterations: config.max_iterations,
            residual: self.get_residual(&state),
        })
    }

    /// Initialize the solution state
    fn initialize(&self, problem: &Self::Problem) -> Result<Self::Solution>;

    /// Get the current residual
    fn get_residual(&self, state: &Self::Solution) -> f64;
}

/// Direct solver trait
pub trait DirectSolver<T: RealField>: Solver<T> {
    /// Assemble the system matrix
    fn assemble(&self, problem: &Self::Problem) -> Result<SystemMatrix<T>>;

    /// Solve the linear system
    fn solve_system(&self, system: &SystemMatrix<T>) -> Result<Self::Solution>;
}

/// System matrix representation
pub struct SystemMatrix<T: RealField> {
    /// Sparse matrix
    pub matrix: nalgebra_sparse::CsrMatrix<T>,
    /// Right-hand side vector
    pub rhs: nalgebra::DVector<T>,
}

/// Solution monitoring trait
pub trait SolutionMonitor<T: RealField> {
    /// Solution type to monitor
    type Solution;

    /// Called at each iteration
    fn on_iteration(&mut self, iteration: usize, solution: &Self::Solution, residual: T);

    /// Called when converged
    fn on_converged(&mut self, iterations: usize, solution: &Self::Solution, residual: T);

    /// Called on failure
    fn on_failed(&mut self, iterations: usize, solution: &Self::Solution, residual: T);
}

/// Default monitor that does nothing
pub struct NullMonitor;

impl<T: RealField> SolutionMonitor<T> for NullMonitor {
    type Solution = ();

    fn on_iteration(&mut self, _: usize, _: &Self::Solution, _: T) {}
    fn on_converged(&mut self, _: usize, _: &Self::Solution, _: T) {}
    fn on_failed(&mut self, _: usize, _: &Self::Solution, _: T) {}
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solver_config() {
        let config = SolverConfig::<f64>::default();
        assert_eq!(config.tolerance, 1e-6);
        assert_eq!(config.max_iterations, 1000);
        assert_eq!(config.relaxation_factor, 1.0);
        assert!(config.parallel);
    }
}