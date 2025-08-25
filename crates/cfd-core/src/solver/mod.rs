//! Solver module providing core solver abstractions
//!
//! This module follows SOLID principles with proper domain separation

pub mod config;
pub mod convergence;
pub mod direct;
pub mod iterative;
pub mod monitor;
pub mod traits;

// Re-export core traits
pub use config::{
    ConvergenceConfig, ExecutionConfig, LinearSolverConfig, NetworkConfig, NetworkSolverConfig,
    NumericalConfig, SolverConfig, SolverConfigBuilder, SolverConfiguration,
};
pub use convergence::{AndCriteria, ConvergenceCriteria, OrCriteria, ToleranceCriteria};
pub use direct::DirectSolver;
pub use iterative::{IterationState, IterativeSolver, SolverIterator};
pub use monitor::{MonitoredIterator, NullMonitor, SolutionMonitor};
pub use traits::{Configurable, Solver, Validatable};
