//! Performance benchmarking framework for CFD operations
//!
//! Provides comprehensive benchmarking capabilities including:
//! - Timing measurements for CFD operations
//! - Memory usage profiling
//! - Parallel scaling analysis
//! - Performance regression detection
//! - Automated benchmarking suites

pub mod analysis;
pub mod memory;
pub mod performance;
pub mod production;
pub mod scaling;
pub mod suite;
pub mod visualization;

pub use analysis::{
    AlertSeverity, PerformanceAnalyzer, PerformanceReport, PerformanceTrend, RegressionAlert,
    RegressionConfig, TrendType,
};
pub use memory::{MemoryProfiler, MemoryStats};
pub use performance::{
    AlgorithmComplexity, CfdPerformanceBenchmarks, PerformanceBenchmark, PerformanceProfile,
    TimingResult,
};
pub use scaling::{ScalingAnalysis, ScalingResult};
pub use suite::{BenchmarkConfig, BenchmarkResult, BenchmarkStatus, BenchmarkSuite};

/// Common benchmarking utilities
pub mod utils {
    use std::time::{Duration, Instant};

    /// High-precision timer for benchmarking
    #[derive(Debug, Clone)]
    pub struct BenchmarkTimer {
        start: Option<Instant>,
        laps: Vec<Duration>,
    }

    impl BenchmarkTimer {
        /// Create a new timer
        pub fn new() -> Self {
            Self {
                start: None,
                laps: Vec::new(),
            }
        }

        /// Start the timer
        pub fn start(&mut self) {
            self.start = Some(Instant::now());
            self.laps.clear();
        }

        /// Record a lap time
        pub fn lap(&mut self) -> Duration {
            if let Some(start) = self.start {
                let elapsed = start.elapsed();
                self.laps.push(elapsed);
                elapsed
            } else {
                Duration::ZERO
            }
        }

        /// Stop the timer and return total elapsed time
        pub fn stop(&mut self) -> Duration {
            if let Some(start) = self.start {
                let total = start.elapsed();
                self.laps.push(total);
                self.start = None;
                total
            } else {
                Duration::ZERO
            }
        }

        /// Get all lap times
        pub fn laps(&self) -> &[Duration] {
            &self.laps
        }

        /// Get the last lap time
        pub fn last_lap(&self) -> Option<Duration> {
            self.laps.last().copied()
        }
    }

    impl Default for BenchmarkTimer {
        fn default() -> Self {
            Self::new()
        }
    }

    /// Statistical analysis for benchmark results using robust statistical measures
    ///
    /// Provides comprehensive statistical characterization of CFD benchmark performance
    /// including central tendency, dispersion, and sample characteristics. All statistics
    /// are computed using numerically stable algorithms suitable for performance analysis.
    ///
    /// # Statistical Properties
    ///
    /// - **Mean**: Arithmetic average of all measurements, represents typical performance
    /// - **Standard Deviation**: Measure of performance variability, lower values indicate more consistent results
    /// - **Min/Max**: Range of observed performance values, useful for outlier detection
    /// - **Median**: Middle value when measurements are sorted, robust to outliers
    /// - **Sample Count**: Number of independent measurements used for statistical analysis
    ///
    /// # Stability Assessment
    ///
    /// The coefficient of variation (CV = std_dev/mean) is used to assess measurement stability.
    /// Lower CV values indicate more reliable and reproducible benchmark results.
    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    pub struct BenchmarkStats {
        /// Arithmetic mean of all benchmark measurements [seconds]
        ///
        /// Represents the expected performance value. Computed as sum(measurements)/count.
        /// For timing benchmarks, this represents the average execution time.
        pub mean: f64,

        /// Standard deviation of benchmark measurements [seconds]
        ///
        /// Measures the dispersion of timing results around the mean. Computed using
        /// sample standard deviation formula: sqrt(sum((x-mean)²)/(n-1)).
        /// Lower values indicate more consistent performance.
        pub std_dev: f64,

        /// Minimum observed measurement value [seconds]
        ///
        /// The fastest execution time recorded. Useful for identifying best-case performance
        /// and detecting measurement outliers.
        pub min: f64,

        /// Maximum observed measurement value [seconds]
        ///
        /// The slowest execution time recorded. Useful for identifying worst-case performance
        /// and detecting measurement outliers or system interference.
        pub max: f64,

        /// Median of benchmark measurements [seconds]
        ///
        /// The middle value when all measurements are sorted. More robust to outliers
        /// than the mean, representing typical performance for skewed distributions.
        /// For even sample counts, computed as average of two middle values.
        pub median: f64,

        /// Number of independent measurements used for statistical analysis
        ///
        /// Total count of benchmark iterations performed. Larger sample sizes provide
        /// more reliable statistical estimates but require longer benchmark execution.
        pub samples: usize,
    }

    impl BenchmarkStats {
        /// Compute statistics from a set of measurements
        pub fn from_measurements(measurements: &[f64]) -> Self {
            if measurements.is_empty() {
                return Self {
                    mean: 0.0,
                    std_dev: 0.0,
                    min: 0.0,
                    max: 0.0,
                    median: 0.0,
                    samples: 0,
                };
            }

            let samples = measurements.len();
            let mean = measurements.iter().sum::<f64>() / samples as f64;

            let variance = measurements.iter().map(|x| (x - mean).powi(2)).sum::<f64>()
                / (samples.saturating_sub(1)) as f64;
            let std_dev = variance.sqrt();

            let mut sorted = measurements.to_vec();
            sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

            let min = *sorted.first().unwrap();
            let max = *sorted.last().unwrap();
            let median = if samples % 2 == 0 {
                (sorted[samples / 2 - 1] + sorted[samples / 2]) / 2.0
            } else {
                sorted[samples / 2]
            };

            Self {
                mean,
                std_dev,
                min,
                max,
                median,
                samples,
            }
        }

        /// Check if measurements are within acceptable variation
        pub fn is_stable(&self, max_cv: f64) -> bool {
            if self.mean == 0.0 {
                return true;
            }
            let cv = self.std_dev / self.mean;
            cv < max_cv
        }

        /// Compute coefficient of variation (CV = std_dev / mean)
        pub fn coefficient_of_variation(&self) -> f64 {
            if self.mean == 0.0 {
                0.0
            } else {
                self.std_dev / self.mean
            }
        }
    }
}

/// Comprehensive benchmark result aggregation and analysis report
///
/// Aggregates all benchmark results from a complete CFD benchmarking session,
/// providing both detailed individual results and high-level summary statistics.
/// Used for performance tracking, regression detection, and optimization validation.
///
/// # Report Components
///
/// - **Individual Results**: Detailed performance data for each benchmark operation
/// - **Summary Statistics**: Aggregated metrics across all benchmarks
/// - **Temporal Context**: Timestamp for tracking performance evolution over time
/// - **Stability Analysis**: Performance consistency assessment across benchmark runs
#[derive(Debug, Clone)]
pub struct BenchmarkReport {
    /// Benchmark suite identifier/name for tracking and comparison
    ///
    /// Human-readable name describing the benchmark suite (e.g., "CFD-2D Navier-Stokes",
    /// "LES Turbulence Models"). Used for report identification and historical comparison.
    pub name: String,

    /// UTC timestamp when the benchmark suite was executed
    ///
    /// Records when the benchmark was performed, enabling temporal analysis of
    /// performance changes, regression detection, and performance trend identification.
    /// Stored in UTC for consistent global comparison.
    pub timestamp: chrono::DateTime<chrono::Utc>,

    /// Collection of individual benchmark results from the suite execution
    ///
    /// Contains detailed performance data for each benchmark operation, including
    /// timing statistics, memory usage, and stability metrics. Each result corresponds
    /// to a specific CFD algorithm or computational kernel being benchmarked.
    pub results: Vec<BenchmarkResult>,

    /// High-level summary statistics across all benchmark results
    ///
    /// Aggregated performance metrics providing an overview of the entire benchmark suite,
    /// including total operations, average throughput, memory usage, and overall stability.
    pub summary: BenchmarkSummary,
}

/// High-level performance summary across an entire benchmark suite
///
/// Provides aggregated performance metrics and stability assessment for a complete
/// CFD benchmarking session. Enables quick evaluation of overall system performance,
/// resource utilization, and result reliability.
///
/// # Key Metrics
///
/// - **Total Operations**: Cumulative count of all computational operations performed
/// - **Total Time**: Wall-clock time for complete benchmark suite execution
/// - **Average Throughput**: Operations per second across all benchmarks
/// - **Memory Usage**: Peak memory consumption during benchmark execution
/// - **Stability Score**: Performance consistency measure (0.0-1.0, higher is better)
///
/// # Performance Assessment
///
/// The stability score combines measurement variability and system consistency
/// to provide a single metric for benchmark reliability assessment.
#[derive(Debug, Clone)]
pub struct BenchmarkSummary {
    /// Total number of computational operations performed across all benchmarks
    ///
    /// Cumulative count of all CFD operations (matrix multiplications, FFTs, etc.)
    /// executed during the benchmark suite. Used for throughput calculations and
    /// computational intensity assessment.
    pub total_operations: usize,

    /// Total wall-clock time for complete benchmark suite execution
    ///
    /// End-to-end execution time from benchmark start to finish, including setup,
    /// computation, and cleanup phases. Used for overall performance assessment
    /// and comparison across different hardware/software configurations.
    pub total_time: std::time::Duration,

    /// Average computational throughput across all benchmark operations [ops/second]
    ///
    /// Computed as total_operations / total_time_seconds. Represents the average
    /// rate of CFD operations performed, providing a standardized performance metric
    /// for comparing different implementations and hardware configurations.
    pub avg_throughput: f64,

    /// Peak memory usage during benchmark execution [bytes]
    ///
    /// Maximum memory consumption observed during the benchmark suite. None if
    /// memory profiling was not enabled. Critical for assessing memory efficiency
    /// and identifying memory bottlenecks in CFD algorithms.
    pub memory_peak: Option<u64>,

    /// Performance stability score (0.0 to 1.0, higher is more stable)
    ///
    /// Composite metric assessing benchmark result consistency and system stability.
    /// Computed from coefficient of variation across all measurements, where 1.0
    /// represents perfectly stable performance and lower values indicate variability.
    /// Used for reliability assessment and outlier detection.
    pub stability_score: f64,
}

impl BenchmarkReport {
    /// Generate a human-readable report
    pub fn format_report(&self) -> String {
        use std::fmt::Write;

        let mut report = String::new();

        writeln!(report, "CFD Performance Benchmark Report").unwrap();
        writeln!(report, "==================================").unwrap();
        writeln!(report, "Benchmark: {}", self.name).unwrap();
        writeln!(
            report,
            "Timestamp: {}",
            self.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
        )
        .unwrap();
        writeln!(
            report,
            "Total Operations: {}",
            self.summary.total_operations
        )
        .unwrap();
        writeln!(
            report,
            "Total Time: {:.3}s",
            self.summary.total_time.as_secs_f64()
        )
        .unwrap();
        writeln!(
            report,
            "Average Throughput: {:.2} ops/s",
            self.summary.avg_throughput
        )
        .unwrap();

        if let Some(memory) = self.summary.memory_peak {
            writeln!(
                report,
                "Peak Memory Usage: {:.2} MB",
                memory as f64 / 1_048_576.0
            )
            .unwrap();
        }

        writeln!(
            report,
            "Stability Score: {:.3}",
            self.summary.stability_score
        )
        .unwrap();
        writeln!(report).unwrap();

        writeln!(report, "Individual Results:").unwrap();
        for result in &self.results {
            if let Some(perf) = &result.performance {
                writeln!(
                    report,
                    "  {}: {:.3}ms ± {:.3}ms ({} samples)",
                    perf.operation_name,
                    perf.stats.mean * 1000.0,
                    perf.stats.std_dev * 1000.0,
                    perf.stats.samples
                )
                .unwrap();
            } else {
                writeln!(report, "  {}: No performance data", result.name).unwrap();
            }
        }

        report
    }
}

#[cfg(test)]
mod tests {
    use super::utils::{BenchmarkStats, BenchmarkTimer};
    use std::time::Duration;

    #[test]
    fn test_benchmark_timer() {
        let mut timer = BenchmarkTimer::new();

        timer.start();
        std::thread::sleep(Duration::from_millis(10));
        let lap1 = timer.lap();

        std::thread::sleep(Duration::from_millis(5));
        let total = timer.stop();

        assert!(lap1 >= Duration::from_millis(10));
        assert!(total >= Duration::from_millis(15));
        assert_eq!(timer.laps().len(), 2);
    }

    #[test]
    fn test_benchmark_stats() {
        let measurements = vec![100.0, 105.0, 95.0, 102.0, 98.0];
        let stats = BenchmarkStats::from_measurements(&measurements);

        assert!((stats.mean - 100.0).abs() < 1.0);
        assert!(stats.std_dev > 0.0 && stats.std_dev < 5.0);
        assert_eq!(stats.samples, 5);
        assert!(stats.is_stable(0.1)); // Should be stable
    }
}
