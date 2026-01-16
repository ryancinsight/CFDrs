//! Performance monitoring and adaptive optimization for CFD operations
//!
//! This module provides runtime performance monitoring to make intelligent
//! decisions about algorithm selection, threshold tuning, and resource allocation.

use std::collections::HashMap;
use std::time::{Duration, Instant};
use std::sync::{Arc, Mutex};

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
    pub fn record_measurement(
        &self,
        operation: &str,
        array_size: usize,
        execution_time_ns: u64,
    ) {
        let mut metrics = self.metrics.lock().unwrap();
        let op_metrics = metrics.entry(operation.to_string()).or_default();
        let size_metrics = op_metrics.entry(array_size).or_default();
        
        // Update running statistics
        let new_count = size_metrics.sample_count + 1;
        let new_avg = (size_metrics.avg_time_ns * size_metrics.sample_count + execution_time_ns) / new_count;
        
        // Simple online variance calculation
        let variance = if new_count > 1 {
            let old_variance = size_metrics.std_dev_ns.pow(2) as f64;
            let diff = execution_time_ns as f64 - size_metrics.avg_time_ns as f64;
            ((new_count - 1) as f64 * old_variance + diff * diff / new_count as f64) / new_count as f64
        } else {
            0.0
        };
        
        size_metrics.avg_time_ns = new_avg;
        size_metrics.std_dev_ns = variance.sqrt() as u64;
        size_metrics.sample_count = new_count;
        size_metrics.throughput_ops_per_sec = 1_000_000_000.0 / new_avg as f64;
        
        // TODO: Implement cache efficiency estimation based on memory access patterns
    }
    
    /// Get optimal threshold for switching between algorithms
    pub fn get_optimal_threshold(&self, operation: &str) -> Option<usize> {
        let metrics = self.metrics.lock().unwrap();
        let op_metrics = metrics.get(operation)?;
        
        // Find crossover point where SIMD becomes beneficial
        let mut sizes: Vec<_> = op_metrics.keys().copied().collect();
        sizes.sort_unstable();
        
        // TODO: Implement more sophisticated threshold detection using regression
        for i in 1..sizes.len() {
            let prev_size = sizes[i - 1];
            let curr_size = sizes[i];
            
            if let (Some(prev_metrics), Some(curr_metrics)) = 
                (op_metrics.get(&prev_size), op_metrics.get(&curr_size)) {
                // Simple heuristic: look for performance improvement
                if curr_metrics.throughput_ops_per_sec > prev_metrics.throughput_ops_per_sec * 1.1 {
                    return Some(curr_size);
                }
            }
        }
        
        None
    }
    
    /// Check if calibration is needed
    pub fn needs_calibration(&self) -> bool {
        let last_calibration = *self.last_calibration.lock().unwrap();
        last_calibration.elapsed() > self.calibration_interval
    }
    
    /// Run performance calibration suite
    /// TODO: Implement comprehensive calibration for all operations
    pub fn run_calibration(&self) -> Result<(), Box<dyn std::error::Error>> {
        println!("Running performance calibration...");
        
        // Mark calibration time
        *self.last_calibration.lock().unwrap() = Instant::now();
        
        // TODO: Implement calibration benchmarks for different operations
        // This should test SIMD vs scalar vs parallel for various array sizes
        
        println!("Performance calibration completed");
        Ok(())
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
            .map(|&threshold| array_size >= threshold)
            .unwrap_or(true) // Default to SIMD if no threshold known
    }
    
    /// Record operation performance for learning
    pub fn record_performance(&mut self, operation: &str, array_size: usize, time_ns: u64) {
        self.monitor.record_measurement(operation, array_size, time_ns);
        
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
