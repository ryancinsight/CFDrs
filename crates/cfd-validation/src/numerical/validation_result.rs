//! Validation result types for numerical algorithms

use nalgebra::{DVector, RealField};

/// Validation result for a numerical algorithm
#[derive(Debug, Clone)]
pub struct ValidationResult<T: RealField + Copy> {
    /// Algorithm name
    pub algorithm_name: String,
    /// Test case name
    pub test_case: String,
    /// Computed solution
    pub computed_solution: DVector<T>,
    /// Analytical solution
    pub analytical_solution: DVector<T>,
    /// Error metrics
    pub error_metrics: super::error_metrics::ErrorMetrics<T>,
    /// Convergence information
    pub convergence_info: ConvergenceInfo<T>,
    /// Literature reference
    pub literature_reference: String,
    /// Test passed
    pub passed: bool,
}

/// Convergence information
#[derive(Debug, Clone)]
pub struct ConvergenceInfo<T: RealField + Copy> {
    /// Number of iterations
    pub iterations: usize,
    /// Final residual norm
    pub final_residual: T,
    /// Convergence rate (if available)
    pub convergence_rate: Option<T>,
}
