use super::summary::{PerformanceMetrics, ValidationSummary};
use cfd_core::error::Result;
use std::collections::HashMap;
use std::time::SystemTime;

/// Main validation report structure
#[derive(Debug, Clone)]
pub struct ValidationReport {
    /// Report generation timestamp
    pub timestamp: SystemTime,
    /// Report title
    pub title: String,
    /// Executive summary
    pub summary: ValidationSummary,
    /// Test results by category
    pub test_results: HashMap<String, TestCategory>,
    /// Performance benchmarks
    pub performance: Vec<PerformanceReport>,
    /// Code quality metrics
    pub code_quality: CodeQualityReport,
    /// Recommendations for improvement
    pub recommendations: Vec<String>,
}

/// Test category results
#[derive(Debug, Clone)]
pub struct TestCategory {
    /// Name of the test category
    pub name: String,
    /// Number of passed tests
    pub passed: usize,
    /// Number of failed tests
    pub failed: usize,
    /// Number of skipped tests
    pub skipped: usize,
    /// Total number of tests
    pub total: usize,
    /// Test coverage percentage
    pub coverage_percentage: f64,
    /// Detailed test results
    pub details: Vec<TestResult>,
}

/// Individual test result
#[derive(Debug, Clone)]
pub struct TestResult {
    /// Name of the test
    pub name: String,
    /// Execution status
    pub status: TestStatus,
    /// Execution duration in milliseconds
    pub duration_ms: f64,
    /// Error message if failed
    pub error_message: Option<String>,
    /// Coverage data if available
    pub coverage_data: Option<String>,
}

/// Test execution status
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TestStatus {
    /// Test passed
    Passed,
    /// Test failed
    Failed,
    /// Test was skipped
    Skipped,
    /// Test timed out
    Timeout,
}

/// Performance benchmark report
#[derive(Debug, Clone)]
pub struct PerformanceReport {
    /// Name of the benchmark
    pub benchmark_name: String,
    /// Performance metrics
    pub metrics: PerformanceMetrics,
    /// Percentage change if regression detected
    pub regression_detected: Option<f64>,
    /// Baseline metrics for comparison
    pub baseline_comparison: Option<PerformanceMetrics>,
}

/// Code quality report
#[derive(Debug, Clone)]
pub struct CodeQualityReport {
    /// Total lines of code
    pub lines_of_code: usize,
    /// Test coverage percentage
    pub test_coverage: f64,
    /// Documentation coverage percentage
    pub documentation_coverage: f64,
    /// Number of clippy warnings
    pub clippy_warnings: usize,
    /// Number of compiler errors
    pub compiler_errors: usize,
    /// Cyclomatic complexity score
    pub cyclomatic_complexity: f64,
    /// Maintainability index
    pub maintainability_index: f64,
}

/// Reporter trait for different output formats
pub trait Reporter {
    /// Generate report string from validation report data
    fn generate_report(&self, report: &ValidationReport) -> Result<String>;
}

impl Default for CodeQualityReport {
    fn default() -> Self {
        Self {
            lines_of_code: 0,
            test_coverage: 0.0,
            documentation_coverage: 0.0,
            clippy_warnings: 0,
            compiler_errors: 0,
            cyclomatic_complexity: 0.0,
            maintainability_index: 0.0,
        }
    }
}

impl ValidationReport {
    /// Generate report in specified format
    pub fn generate<T: Reporter>(&self, reporter: &T) -> Result<String> {
        reporter.generate_report(self)
    }

    /// Calculate overall health score (0.0 to 1.0)
    #[must_use]
    pub fn health_score(&self) -> f64 {
        let test_score = if self.summary.total_tests > 0 {
            self.summary.passed_tests as f64 / self.summary.total_tests as f64
        } else {
            0.0
        };

        let coverage_score = self.code_quality.test_coverage / 100.0;
        let quality_score = 1.0 - (self.code_quality.clippy_warnings as f64 * 0.01).min(1.0);

        0.4 * test_score + 0.4 * coverage_score + 0.2 * quality_score
    }

    /// Get critical issues requiring attention
    #[must_use]
    pub fn critical_issues(&self) -> Vec<String> {
        let mut issues = Vec::new();

        if self.summary.failed_tests > 0 {
            issues.push(format!(
                "{} test failures detected",
                self.summary.failed_tests
            ));
        }

        if self.code_quality.test_coverage < 80.0 {
            issues.push(format!(
                "Test coverage below 80%: {:.1}%",
                self.code_quality.test_coverage
            ));
        }

        if self.code_quality.clippy_warnings > 10 {
            issues.push(format!(
                "High number of clippy warnings: {}",
                self.code_quality.clippy_warnings
            ));
        }

        let regressions: Vec<_> = self
            .performance
            .iter()
            .filter(|p| p.regression_detected.is_some())
            .collect();

        if !regressions.is_empty() {
            issues.push(format!(
                "{} performance regressions detected",
                regressions.len()
            ));
        }

        issues
    }
}
