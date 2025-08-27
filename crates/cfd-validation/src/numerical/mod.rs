//! Numerical algorithm validation module
//!
//! This module provides comprehensive validation tests for numerical algorithms
//! implemented in the CFD suite, comparing against known analytical solutions
//! and published benchmark results.

pub mod error_metrics;
pub mod linear_solver;
pub mod test_cases;
pub mod validation_result;

// Re-export main types
pub use error_metrics::{compute_error_metrics, ErrorMetrics};
pub use linear_solver::LinearSolverValidator;
pub use validation_result::{ConvergenceInfo, ValidationResult};
