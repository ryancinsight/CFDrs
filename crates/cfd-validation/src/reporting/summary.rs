//! Validation summary structures

use std::time::Duration;

/// Validation test summary
#[derive(Debug, Clone)]
pub struct ValidationSummary {
    /// Total number of tests executed
    pub total_tests: usize,
    /// Number of passed tests
    pub passed_tests: usize,
    /// Number of failed tests
    pub failed_tests: usize,
    /// Number of skipped tests
    pub skipped_tests: usize,
    /// Total execution duration
    pub total_duration: Duration,
    /// Test coverage percentage
    pub coverage_percentage: f64,
}

impl Default for ValidationSummary {
    fn default() -> Self {
        Self {
            total_tests: 0,
            passed_tests: 0,
            failed_tests: 0,
            skipped_tests: 0,
            total_duration: Duration::default(),
            coverage_percentage: 0.0,
        }
    }
}

impl ValidationSummary {
    /// Calculate the percentage of tests that passed
    pub fn pass_rate(&self) -> f64 {
        if self.total_tests == 0 {
            0.0
        } else {
            self.passed_tests as f64 / self.total_tests as f64
        }
    }

    /// Calculate the percentage of tests that failed
    pub fn failure_rate(&self) -> f64 {
        if self.total_tests == 0 {
            0.0
        } else {
            self.failed_tests as f64 / self.total_tests as f64
        }
    }
}

/// Performance metrics structure
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceMetrics {
    /// Mean execution time
    pub mean: f64,
    /// Standard deviation of execution time
    pub std_dev: f64,
    /// Minimum execution time
    pub min: f64,
    /// Maximum execution time
    pub max: f64,
    /// Median execution time
    pub median: f64,
    /// Number of samples collected
    pub samples: usize,
}

impl PerformanceMetrics {
    /// Check if performance is stable (std_dev / mean < threshold)
    pub fn is_stable(&self, threshold: f64) -> bool {
        self.std_dev / self.mean < threshold
    }

    /// Calculate coefficient of variation (std_dev / mean)
    pub fn coefficient_of_variation(&self) -> f64 {
        if self.mean == 0.0 {
            0.0
        } else {
            self.std_dev / self.mean
        }
    }
}

impl Default for PerformanceMetrics {
    fn default() -> Self {
        Self {
            mean: 0.0,
            std_dev: 0.0,
            min: 0.0,
            max: 0.0,
            median: 0.0,
            samples: 0,
        }
    }
}
