//! Performance monitoring and adaptive optimization for CFD operations
//!
//! This module provides runtime performance monitoring to make intelligent
//! decisions about algorithm selection, threshold tuning, and resource allocation.

use std::collections::HashMap;
use std::hint::black_box;
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

/// Performance metrics for different operation types
#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    /// Average execution time in nanoseconds
    pub avg_time_ns: u64,
    /// Standard deviation of execution times
    pub std_dev_ns: u64,
    /// Number of samples collected
    pub sample_count: u64,
    /// Throughput in operations per second
    pub throughput_ops_per_sec: f64,
    /// Cache efficiency estimate (0.0 to 1.0)
    pub cache_efficiency: f64,
}

impl Default for PerformanceMetrics {
    fn default() -> Self {
        Self {
            avg_time_ns: 0,
            std_dev_ns: 0,
            sample_count: 0,
            throughput_ops_per_sec: 0.0,
            cache_efficiency: 0.0,
        }
    }
}

/// Adaptive performance monitor that learns optimal thresholds
pub struct AdaptivePerformanceMonitor {
    /// Performance history for different operations and sizes
    metrics: Arc<Mutex<HashMap<String, HashMap<usize, PerformanceMetrics>>>>,
    /// Last calibration time
    last_calibration: Arc<Mutex<Instant>>,
    /// Calibration interval
    calibration_interval: Duration,
}

impl AdaptivePerformanceMonitor {
    /// Create new performance monitor
    pub fn new() -> Self {
        Self {
            metrics: Arc::new(Mutex::new(HashMap::new())),
            last_calibration: Arc::new(Mutex::new(Instant::now())),
            calibration_interval: Duration::from_secs(60), // Calibrate every minute
        }
    }

    /// Record performance measurement for an operation
    pub fn record_measurement(&self, operation: &str, array_size: usize, execution_time_ns: u64) {
        let mut metrics = self.metrics.lock().unwrap();
        let op_metrics = metrics.entry(operation.to_string()).or_default();
        let (old_avg, old_std, old_count) = {
            let size_metrics = op_metrics.entry(array_size).or_default();
            (
                size_metrics.avg_time_ns,
                size_metrics.std_dev_ns,
                size_metrics.sample_count,
            )
        };

        let new_count = old_count + 1;
        let new_avg = (old_avg * old_count + execution_time_ns) / new_count;

        let variance = if new_count > 1 {
            let old_variance = old_std.pow(2) as f64;
            let diff = execution_time_ns as f64 - old_avg as f64;
            ((new_count - 1) as f64 * old_variance + diff * diff / new_count as f64)
                / new_count as f64
        } else {
            0.0
        };

        {
            let size_metrics = op_metrics.entry(array_size).or_default();
            size_metrics.avg_time_ns = new_avg;
            size_metrics.std_dev_ns = variance.sqrt() as u64;
            size_metrics.sample_count = new_count;
            size_metrics.throughput_ops_per_sec = 1_000_000_000.0 / new_avg as f64;
        }

        let max_throughput = op_metrics
            .values()
            .map(|metrics| metrics.throughput_ops_per_sec)
            .fold(0.0_f64, f64::max);
        let cv = if new_avg > 0 {
            variance.sqrt() / new_avg as f64
        } else {
            0.0
        };
        let stability = 1.0 / (1.0 + cv);
        let throughput_ops_per_sec = 1_000_000_000.0 / new_avg as f64;
        let relative_throughput = if max_throughput > 0.0 {
            throughput_ops_per_sec / max_throughput
        } else {
            0.0
        };
        let cache_efficiency = (relative_throughput * stability).clamp(0.0, 1.0);
        let size_metrics = op_metrics.entry(array_size).or_default();
        size_metrics.cache_efficiency = cache_efficiency;
    }

    /// Get optimal threshold for switching between algorithms
    pub fn get_optimal_threshold(&self, operation: &str) -> Option<usize> {
        let metrics = self.metrics.lock().unwrap();
        let op_metrics = metrics.get(operation)?;

        // Find crossover point where SIMD becomes beneficial
        let mut sizes: Vec<_> = op_metrics.keys().copied().collect();
        sizes.sort_unstable();

        let mut x = Vec::with_capacity(sizes.len());
        let mut y = Vec::with_capacity(sizes.len());
        for size in &sizes {
            let throughput = op_metrics.get(size)?.throughput_ops_per_sec;
            x.push(*size as f64);
            y.push(throughput);
        }

        if sizes.len() < 4 {
            return None;
        }

        let mut best_split = None;
        let mut best_sse = f64::INFINITY;
        let mut best_slopes = (0.0, 0.0);

        for split in 2..sizes.len() - 1 {
            let (slope_left, intercept_left, sse_left) =
                Self::linear_regression(&x[..split], &y[..split])?;
            let (slope_right, intercept_right, sse_right) =
                Self::linear_regression(&x[split..], &y[split..])?;
            let total_sse = sse_left + sse_right;
            if total_sse < best_sse {
                best_sse = total_sse;
                best_split = Some(split);
                best_slopes = (slope_left, slope_right);
                let _ = intercept_left;
                let _ = intercept_right;
            }
        }

        let split = best_split?;
        if best_slopes.1 > best_slopes.0 * 1.05 {
            Some(sizes[split])
        } else {
            None
        }
    }

    /// Check if calibration is needed
    pub fn needs_calibration(&self) -> bool {
        let last_calibration = *self.last_calibration.lock().unwrap();
        last_calibration.elapsed() > self.calibration_interval
    }

    /// Run performance calibration suite
    pub fn run_calibration(&self) -> Result<(), Box<dyn std::error::Error>> {
        println!("Running performance calibration...");
        let sizes = [64usize, 128, 256, 512, 1024, 2048, 4096, 8192];

        for size in sizes {
            let iterations = (100_000usize / size.max(1)).max(10);
            let a = vec![1.0_f64; size];
            let b = vec![2.0_f64; size];
            let mut c = vec![0.0_f64; size];
            let alpha = 0.5_f64;

            let start = Instant::now();
            for _ in 0..iterations {
                for i in 0..size {
                    c[i] = black_box(a[i] + b[i]);
                }
            }
            let elapsed = start.elapsed();
            let avg = elapsed.as_nanos() as u64 / iterations as u64;
            self.record_measurement("vector_add", size, avg);

            let start = Instant::now();
            for _ in 0..iterations {
                for i in 0..size {
                    c[i] = black_box(a[i] * b[i]);
                }
            }
            let elapsed = start.elapsed();
            let avg = elapsed.as_nanos() as u64 / iterations as u64;
            self.record_measurement("vector_mul", size, avg);

            let start = Instant::now();
            for _ in 0..iterations {
                let mut sum = 0.0_f64;
                for i in 0..size {
                    sum += black_box(a[i] * b[i]);
                }
                black_box(sum);
            }
            let elapsed = start.elapsed();
            let avg = elapsed.as_nanos() as u64 / iterations as u64;
            self.record_measurement("dot", size, avg);

            let start = Instant::now();
            for _ in 0..iterations {
                for i in 0..size {
                    c[i] = black_box(alpha * a[i] + b[i]);
                }
            }
            let elapsed = start.elapsed();
            let avg = elapsed.as_nanos() as u64 / iterations as u64;
            self.record_measurement("axpy", size, avg);

            let nnz_per_row = 3usize.min(size.max(1));
            let mut row_offsets = Vec::with_capacity(size + 1);
            let mut col_indices = Vec::with_capacity(size * nnz_per_row);
            let mut values = Vec::with_capacity(size * nnz_per_row);
            row_offsets.push(0u32);
            for i in 0..size {
                let mut cols = Vec::new();
                if i > 0 {
                    cols.push(i - 1);
                }
                cols.push(i);
                if i + 1 < size {
                    cols.push(i + 1);
                }
                cols.sort_unstable();
                for col in cols {
                    col_indices.push(col as u32);
                    values.push(1.0_f64 / (col as f64 + 1.0));
                }
                row_offsets.push(col_indices.len() as u32);
            }

            let start = Instant::now();
            for _ in 0..iterations {
                for i in 0..size {
                    let start_idx = row_offsets[i] as usize;
                    let end_idx = row_offsets[i + 1] as usize;
                    let mut sum = 0.0_f64;
                    for idx in start_idx..end_idx {
                        sum += values[idx] * b[col_indices[idx] as usize];
                    }
                    c[i] = black_box(sum);
                }
            }
            let elapsed = start.elapsed();
            let avg = elapsed.as_nanos() as u64 / iterations as u64;
            self.record_measurement("spmv", size, avg);
        }

        *self.last_calibration.lock().unwrap() = Instant::now();

        println!("Performance calibration completed");
        Ok(())
    }

    fn linear_regression(x: &[f64], y: &[f64]) -> Option<(f64, f64, f64)> {
        let n = x.len();
        if n < 2 {
            return None;
        }
        let n_f = n as f64;
        let sum_x = x.iter().sum::<f64>();
        let sum_y = y.iter().sum::<f64>();
        let sum_x2 = x.iter().map(|v| v * v).sum::<f64>();
        let sum_xy = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum::<f64>();
        let denom = n_f * sum_x2 - sum_x * sum_x;
        if denom == 0.0 {
            return None;
        }
        let slope = (n_f * sum_xy - sum_x * sum_y) / denom;
        let intercept = (sum_y - slope * sum_x) / n_f;
        let mut sse = 0.0_f64;
        for (xi, yi) in x.iter().zip(y.iter()) {
            let pred = slope * xi + intercept;
            let resid = yi - pred;
            sse += resid * resid;
        }
        Some((slope, intercept, sse))
    }

    /// Get performance report
    pub fn get_performance_report(&self) -> HashMap<String, HashMap<usize, PerformanceMetrics>> {
        self.metrics.lock().unwrap().clone()
    }
}

impl Default for AdaptivePerformanceMonitor {
    fn default() -> Self {
        Self::new()
    }
}

/// Performance-aware operation selector
pub struct PerformanceAwareSelector {
    monitor: AdaptivePerformanceMonitor,
    default_thresholds: HashMap<String, usize>,
}

impl PerformanceAwareSelector {
    /// Create new selector with default thresholds
    pub fn new() -> Self {
        let mut default_thresholds = HashMap::new();
        default_thresholds.insert("vector_add".to_string(), 500);
        default_thresholds.insert("vector_mul".to_string(), 500);
        default_thresholds.insert("spmv".to_string(), 1000);

        Self {
            monitor: AdaptivePerformanceMonitor::new(),
            default_thresholds,
        }
    }

    /// Decide whether to use SIMD for given operation and size
    pub fn should_use_simd(&mut self, operation: &str, array_size: usize) -> bool {
        // Check if we have learned a better threshold
        if let Some(optimal_threshold) = self.monitor.get_optimal_threshold(operation) {
            return array_size >= optimal_threshold;
        }

        // Fall back to default threshold
        self.default_thresholds
            .get(operation)
            .is_none_or(|&threshold| array_size >= threshold) // Default to SIMD if no threshold known
    }

    /// Record operation performance for learning
    pub fn record_performance(&mut self, operation: &str, array_size: usize, time_ns: u64) {
        self.monitor
            .record_measurement(operation, array_size, time_ns);

        // Run calibration if needed
        if self.monitor.needs_calibration() {
            let _ = self.monitor.run_calibration();
        }
    }

    /// Get current performance monitor
    pub fn monitor(&self) -> &AdaptivePerformanceMonitor {
        &self.monitor
    }
}

impl Default for PerformanceAwareSelector {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_performance_metrics_recording() {
        let monitor = AdaptivePerformanceMonitor::new();

        // Record some measurements
        monitor.record_measurement("test_op", 100, 1000);
        monitor.record_measurement("test_op", 100, 1100);
        monitor.record_measurement("test_op", 100, 900);

        let report = monitor.get_performance_report();
        let metrics = &report["test_op"][&100];

        assert_eq!(metrics.sample_count, 3);
        assert_eq!(metrics.avg_time_ns, 1000); // Average of 1000, 1100, 900
        assert!(metrics.throughput_ops_per_sec > 0.0);
    }

    #[test]
    fn test_performance_aware_selector() {
        let mut selector = PerformanceAwareSelector::new();

        // Should use default threshold initially
        assert!(selector.should_use_simd("vector_add", 1000));
        assert!(!selector.should_use_simd("vector_add", 100));

        // Record performance data
        selector.record_performance("vector_add", 100, 500);
        selector.record_performance("vector_add", 1000, 100);

        // Decision should be based on recorded performance
        assert!(selector.should_use_simd("vector_add", 1000));
    }
}
