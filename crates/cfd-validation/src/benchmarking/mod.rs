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

pub use analysis::{AlertSeverity, PerformanceAnalyzer, PerformanceReport, PerformanceTrend, RegressionAlert, RegressionConfig, TrendType};
pub use memory::{MemoryProfiler, MemoryStats};
pub use performance::{PerformanceBenchmark, TimingResult};
pub use scaling::{ScalingAnalysis, ScalingResult};
pub use suite::{BenchmarkConfig, BenchmarkSuite, BenchmarkResult, BenchmarkStatus};

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

    /// Statistical analysis for benchmark results
    #[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
    pub struct BenchmarkStats {
        pub mean: f64,
        pub std_dev: f64,
        pub min: f64,
        pub max: f64,
        pub median: f64,
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

            let variance = measurements.iter()
                .map(|x| (x - mean).powi(2))
                .sum::<f64>() / (samples.saturating_sub(1)) as f64;
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

/// Benchmark result aggregation and analysis
#[derive(Debug, Clone)]
pub struct BenchmarkReport {
    pub name: String,
    pub timestamp: chrono::DateTime<chrono::Utc>,
    pub results: Vec<BenchmarkResult>,
    pub summary: BenchmarkSummary,
}

#[derive(Debug, Clone)]
pub struct BenchmarkSummary {
    pub total_operations: usize,
    pub total_time: std::time::Duration,
    pub avg_throughput: f64,
    pub memory_peak: Option<u64>,
    pub stability_score: f64, // 0.0 to 1.0, higher is more stable
}

impl BenchmarkReport {
    /// Generate a human-readable report
    pub fn format_report(&self) -> String {
        use std::fmt::Write;

        let mut report = String::new();

        writeln!(report, "CFD Performance Benchmark Report").unwrap();
        writeln!(report, "==================================").unwrap();
        writeln!(report, "Benchmark: {}", self.name).unwrap();
        writeln!(report, "Timestamp: {}", self.timestamp.format("%Y-%m-%d %H:%M:%S UTC")).unwrap();
        writeln!(report, "Total Operations: {}", self.summary.total_operations).unwrap();
        writeln!(report, "Total Time: {:.3}s", self.summary.total_time.as_secs_f64()).unwrap();
        writeln!(report, "Average Throughput: {:.2} ops/s", self.summary.avg_throughput).unwrap();

        if let Some(memory) = self.summary.memory_peak {
            writeln!(report, "Peak Memory Usage: {:.2} MB", memory as f64 / 1_048_576.0).unwrap();
        }

        writeln!(report, "Stability Score: {:.3}", self.summary.stability_score).unwrap();
        writeln!(report).unwrap();

        writeln!(report, "Individual Results:").unwrap();
        for result in &self.results {
            if let Some(perf) = &result.performance {
                writeln!(report, "  {}: {:.3}ms ± {:.3}ms ({} samples)",
                        perf.operation_name,
                        perf.stats.mean * 1000.0,
                        perf.stats.std_dev * 1000.0,
                        perf.stats.samples).unwrap();
            } else {
                writeln!(report, "  {}: No performance data", result.name).unwrap();
            }
        }

        report
    }
}

#[cfg(test)]
mod tests {
    use super::utils::{BenchmarkTimer, BenchmarkStats};
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
