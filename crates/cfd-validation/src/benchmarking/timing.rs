//! Timing utilities for CFD performance benchmarks
//!
//! Domain: CFD Performance Benchmarking
//! Bounded Context: Timing and Profiling
//! Architectural Pattern: Utility

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
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BenchmarkStats {
    /// Arithmetic mean of all benchmark measurements [seconds]
    pub mean: f64,

    /// Standard deviation of benchmark measurements [seconds]
    pub std_dev: f64,

    /// Minimum observed measurement value [seconds]
    pub min: f64,

    /// Maximum observed measurement value [seconds]
    pub max: f64,

    /// Median of benchmark measurements [seconds]
    pub median: f64,

    /// Number of independent measurements used for statistical analysis
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
        let median = if samples.is_multiple_of(2) {
            f64::midpoint(sorted[samples / 2 - 1], sorted[samples / 2])
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

/// Utility for tracking the specific elapsed time of a scoped block
pub struct ProfileSpan {
    name: String,
    start: Instant,
}

impl ProfileSpan {
    /// Start a new profile span
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            start: Instant::now(),
        }
    }
}

impl Drop for ProfileSpan {
    fn drop(&mut self) {
        let elapsed = self.start.elapsed();
        tracing::trace!("ProfileSpan [{}]: {:?}", self.name, elapsed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
