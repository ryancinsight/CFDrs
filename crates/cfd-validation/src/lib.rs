//! CFD validation library
//! Validation and benchmarking tools for CFD simulations.

#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::must_use_candidate)]

pub mod analytical;
pub mod convergence;  // Now a modular directory structure
pub mod error_metrics;
pub mod numerical_validation;
pub mod time_integration_validation;

// Re-export commonly used types
pub use convergence::{
    ConvergenceStudy, 
    RichardsonExtrapolation,
    ConvergenceAnalysis,
    ConvergenceOrder,
    ConvergenceStatus,
    ConvergenceCriterion,
    GridConvergenceIndex,
};
pub use error_metrics::{ErrorMetric, ErrorAnalysis};
pub use numerical_validation::LinearSolverValidator;