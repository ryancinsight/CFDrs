//! Comprehensive benchmarking suite for CFD operations
//!
//! Orchestrates performance, memory, and scaling analysis into
//! automated benchmark suites with reporting and regression detection.
//!
//! Domain: CFD Performance Benchmarking
//! Bounded Context: Benchmark Suite Orchestration
//! Architectural Pattern: Builder + Facade

use super::{memory::MemoryStats, performance::TimingResult, scaling::ScalingResult};
use cfd_core::error::{Error, Result};
use std::collections::HashMap;
use std::time::Duration;

/// Benchmark configuration with architectural purity
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BenchmarkConfig {
    /// Number of iterations for each benchmark
    pub iterations: usize,
    /// Enable memory profiling
    pub enable_memory: bool,
    /// Enable scaling analysis
    pub enable_scaling: bool,
    /// Enable detailed reporting
    pub detailed_reporting: bool,
    /// Performance regression threshold (percentage)
    pub regression_threshold: f64,
    /// Baseline performance data (for regression detection)
    pub baseline_data: Option<HashMap<String, f64>>,
    /// Problem sizes to benchmark
    pub problem_sizes: Vec<usize>,
    /// Output directory for reports
    pub output_dir: Option<String>,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            iterations: 5,
            enable_memory: true,
            enable_scaling: true,
            detailed_reporting: true,
            regression_threshold: 5.0, // 5% regression threshold
            baseline_data: None,
            problem_sizes: vec![1000, 10000, 100000],
            output_dir: Some("performance_results".to_string()),
        }
    }
}

/// Individual benchmark result with comprehensive metrics
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BenchmarkResult {
    /// Benchmark name (operation identifier)
    pub name: String,
    /// Problem size used for this benchmark
    pub problem_size: usize,
    /// Performance timing results
    pub performance: Option<TimingResult>,
    /// Memory usage statistics
    pub memory: Option<MemoryStats>,
    /// Scaling analysis results
    pub scaling: Option<ScalingResult>,
    /// Execution status
    pub status: BenchmarkStatus,
    /// Performance regression detected
    pub regression_detected: Option<f64>, // Percentage change from baseline
    /// Execution duration
    pub duration: Duration,
    /// Timestamp of execution
    pub timestamp: chrono::DateTime<chrono::Utc>,
}

impl BenchmarkResult {
    /// Create new benchmark result with architectural purity
    pub fn new(name: String, problem_size: usize) -> Self {
        Self {
            name,
            problem_size,
            performance: None,
            memory: None,
            scaling: None,
            status: BenchmarkStatus::Skipped,
            regression_detected: None,
            duration: Duration::ZERO,
            timestamp: chrono::Utc::now(),
        }
    }

    /// Set performance metrics
    pub fn with_performance(mut self, timing: TimingResult) -> Self {
        self.performance = Some(timing);
        self
    }

    /// Set memory statistics
    pub fn with_memory(mut self, memory: MemoryStats) -> Self {
        self.memory = Some(memory);
        self
    }

    /// Set scaling results
    pub fn with_scaling(mut self, scaling: ScalingResult) -> Self {
        self.scaling = Some(scaling);
        self
    }

    /// Set execution status
    pub fn with_status(mut self, status: BenchmarkStatus) -> Self {
        self.status = status;
        self
    }

    /// Set regression detection
    pub fn with_regression(mut self, regression: Option<f64>) -> Self {
        self.regression_detected = regression;
        self
    }

    /// Set execution duration
    pub fn with_duration(mut self, duration: Duration) -> Self {
        self.duration = duration;
        self
    }
}

/// Benchmark execution status with architectural enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum BenchmarkStatus {
    /// Benchmark completed successfully
    Passed,
    /// Benchmark failed to execute
    Failed,
    /// Performance regression detected
    Regression,
    /// Benchmark was skipped
    Skipped,
    /// Benchmark is currently running
    Running,
}

/// Comprehensive benchmark suite orchestrator
pub struct BenchmarkSuite {
    config: BenchmarkConfig,
    results: Vec<BenchmarkResult>,
    start_time: std::time::Instant,
}

impl BenchmarkSuite {
    /// Create new benchmark suite with configuration
    pub fn new(config: BenchmarkConfig) -> Self {
        Self {
            config,
            results: Vec::new(),
            start_time: std::time::Instant::now(),
        }
    }

    /// Create with default configuration
    pub fn default() -> Self {
        Self::new(BenchmarkConfig::default())
    }

    /// Add benchmark result to suite
    pub fn add_result(&mut self, result: BenchmarkResult) {
        self.results.push(result);
    }

    /// Get all results
    pub fn results(&self) -> &[BenchmarkResult] {
        &self.results
    }

    /// Get mutable results for modification
    pub fn results_mut(&mut self) -> &mut Vec<BenchmarkResult> {
        &mut self.results
    }

    /// Get configuration
    pub fn config(&self) -> &BenchmarkConfig {
        &self.config
    }

    /// Calculate suite statistics
    pub fn statistics(&self) -> SuiteStatistics {
        let total = self.results.len();
        let passed = self.results.iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Passed))
            .count();
        let failed = self.results.iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Failed))
            .count();
        let regressions = self.results.iter()
            .filter(|r| matches!(r.status, BenchmarkStatus::Regression) || r.regression_detected.is_some())
            .count();

        let total_duration = self.results.iter()
            .map(|r| r.duration)
            .sum();

        SuiteStatistics {
            total_benchmarks: total,
            passed,
            failed,
            regressions,
            total_duration,
            average_duration: if total > 0 { total_duration / total as u32 } else { Duration::ZERO },
        }
    }

    /// Generate comprehensive report
    pub fn generate_report(&self) -> Result<String> {
        let stats = self.statistics();
        let mut report = format!(
            "CFD Benchmark Suite Report\n{}\n\n",
            "=".repeat(50)
        );

        report.push_str(&format!("Configuration:\n"));
        report.push_str(&format!("  Iterations: {}\n", self.config.iterations));
        report.push_str(&format!("  Memory Profiling: {}\n", self.config.enable_memory));
        report.push_str(&format!("  Scaling Analysis: {}\n", self.config.enable_scaling));
        report.push_str(&format!("  Problem Sizes: {:?}\n\n", self.config.problem_sizes));

        report.push_str(&format!("Suite Statistics:\n"));
        report.push_str(&format!("  Total Benchmarks: {}\n", stats.total_benchmarks));
        report.push_str(&format!("  Passed: {} ({:.1}%)\n",
            stats.passed,
            (stats.passed as f64 / stats.total_benchmarks as f64) * 100.0
        ));
        report.push_str(&format!("  Failed: {}\n", stats.failed));
        report.push_str(&format!("  Regressions: {}\n", stats.regressions));
        report.push_str(&format!("  Total Duration: {:.3}s\n", stats.total_duration.as_secs_f64()));
        report.push_str(&format!("  Average Duration: {:.3}s\n\n", stats.average_duration.as_secs_f64()));

        report.push_str("Benchmark Results:\n");
        for result in &self.results {
            let status_icon = match result.status {
                BenchmarkStatus::Passed => "✅",
                BenchmarkStatus::Failed => "❌",
                BenchmarkStatus::Regression => "📉",
                BenchmarkStatus::Skipped => "⏭️",
                BenchmarkStatus::Running => "⏳",
            };

            report.push_str(&format!("  {} {} (size: {}): {:.3}ms",
                status_icon,
                result.name,
                result.problem_size,
                result.duration.as_secs_f64() * 1000.0
            ));

            if let Some(regression) = result.regression_detected {
                report.push_str(&format!(" [Regression: {:.1}%]", regression));
            }
            report.push_str("\n");
        }

        Ok(report)
    }
}

/// Suite-level statistics
#[derive(Debug, Clone)]
pub struct SuiteStatistics {
    pub total_benchmarks: usize,
    pub passed: usize,
    pub failed: usize,
    pub regressions: usize,
    pub total_duration: Duration,
    pub average_duration: Duration,
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_benchmark_result_builder() {
        let result = BenchmarkResult::new("test_op".to_string(), 1000)
            .with_status(BenchmarkStatus::Passed)
            .with_duration(Duration::from_millis(150));

        assert_eq!(result.name, "test_op");
        assert_eq!(result.problem_size, 1000);
        assert_eq!(result.status, BenchmarkStatus::Passed);
        assert_eq!(result.duration, Duration::from_millis(150));
    }

    #[test]
    fn test_suite_statistics() {
        let mut suite = BenchmarkSuite::default();

        // Add test results
        suite.add_result(BenchmarkResult::new("op1".to_string(), 100)
            .with_status(BenchmarkStatus::Passed)
            .with_duration(Duration::from_millis(100)));

        suite.add_result(BenchmarkResult::new("op2".to_string(), 100)
            .with_status(BenchmarkStatus::Failed)
            .with_duration(Duration::from_millis(50)));

        let stats = suite.statistics();
        assert_eq!(stats.total_benchmarks, 2);
        assert_eq!(stats.passed, 1);
        assert_eq!(stats.failed, 1);
        assert_eq!(stats.total_duration, Duration::from_millis(150));
    }

    #[test]
    fn test_config_default() {
        let config = BenchmarkConfig::default();
        assert_eq!(config.iterations, 5);
        assert!(config.enable_memory);
        assert!(config.enable_scaling);
        assert_eq!(config.problem_sizes, vec![1000, 10000, 100000]);
    }
}
