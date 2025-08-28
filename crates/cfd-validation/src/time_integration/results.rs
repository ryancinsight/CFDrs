//! Time integration validation results.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Result of a time integration validation test
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TimeIntegrationResult<T: RealField + Copy> {
    /// Name of the integration method
    pub method_name: String,
    /// Name of the test problem
    pub test_problem: String,
    /// Time step size used
    pub dt: T,
    /// Final time
    pub final_time: T,
    /// Number of time steps
    pub n_steps: usize,
    /// Global error at final time
    pub global_error: T,
    /// Whether the test passed (error below threshold)
    pub passed: bool,
    /// Expected order of accuracy
    pub expected_order: usize,
    /// Computed order of accuracy (if available)
    pub computed_order: Option<T>,
}

impl<T: RealField + Copy> TimeIntegrationResult<T> {
    /// Create a result
    pub fn create(
        method_name: String,
        test_problem: String,
        dt: T,
        final_time: T,
        n_steps: usize,
        global_error: T,
        passed: bool,
        expected_order: usize,
    ) -> Self {
        Self {
            method_name,
            test_problem,
            dt,
            final_time,
            n_steps,
            global_error,
            passed,
            expected_order,
            computed_order: None,
        }
    }

    /// Set computed order of accuracy
    pub fn with_computed_order(mut self, order: T) -> Self {
        self.computed_order = Some(order);
        self
    }
}
