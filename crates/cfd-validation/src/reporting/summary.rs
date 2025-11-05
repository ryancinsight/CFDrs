//! Validation summary structures

use std::time::Duration;

/// Validation test summary
#[derive(Debug, Clone)]
pub struct ValidationSummary {
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub skipped_tests: usize,
    pub total_duration: Duration,
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
    pub fn pass_rate(&self) -> f64 {
        if self.total_tests == 0 {
            0.0
        } else {
            self.passed_tests as f64 / self.total_tests as f64
        }
    }

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
    pub mean: f64,
    pub std_dev: f64,
    pub min: f64,
    pub max: f64,
    pub median: f64,
    pub samples: usize,
}

impl PerformanceMetrics {
    pub fn is_stable(&self, threshold: f64) -> bool {
        self.std_dev / self.mean < threshold
    }

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
