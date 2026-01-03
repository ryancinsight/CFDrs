//! Solver module providing core solver abstractions
//!
//! This module follows SOLID principles with proper domain separation

pub mod config;
pub mod convergence;
pub mod direct;
pub mod iterative;
pub mod monitor;
pub mod traits;

// Curated Re-exports: Only expose the most important, high-level types.
pub use self::config::{SolverConfig, SolverConfiguration};
pub use self::traits::{Configurable, Solver, Validatable};

/// A trait for solvers that compute a solution in a single step (e.g., via matrix inversion).
pub use self::direct::DirectSolver;
/// A trait for solvers that approach a solution through successive approximations.
pub use self::iterative::IterativeSolver;

// Public Sub-modules: Expose specialized functionality in organized namespaces.

/// Solver configuration types.
pub mod configuration {
    pub use super::config::{
        ConvergenceConfig, ExecutionConfig, LinearSolverConfig, NetworkConfig, NetworkSolverConfig,
        NumericalConfig, SolverConfigBuilder,
    };
}

/// Convergence criteria for iterative solvers.
pub mod convergence_criteria {
    pub use super::convergence::{AndCriteria, ConvergenceCriteria, OrCriteria, ToleranceCriteria};
}

/// Solution monitoring tools.
pub mod monitoring {
    pub use super::monitor::{MonitoredIterator, NullMonitor, SolutionMonitor};
}

/// Iteration state management.
pub mod iteration {
    pub use super::iterative::{IterationState, SolverIterator};
}
