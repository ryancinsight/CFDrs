//! Solver module providing core solver abstractions
//!
//! This module follows SOLID principles with proper domain separation

pub mod traits;
pub mod config;
pub mod iterative;
pub mod direct;
pub mod convergence;
pub mod monitor;

// Re-export core traits
pub use traits::{Solver, Configurable, Validatable};
pub use config::{
    SolverConfiguration, SolverConfig, SolverConfigBuilder,
    LinearSolverConfig, NetworkSolverConfig, NetworkConfig,
    ConvergenceConfig, ExecutionConfig, NumericalConfig
};
pub use iterative::{IterativeSolver, IterationState, SolverIterator};
pub use direct::DirectSolver;
pub use convergence::{ConvergenceCriteria, ToleranceCriteria, AndCriteria, OrCriteria};
pub use monitor::{SolutionMonitor, MonitoredIterator, NullMonitor};