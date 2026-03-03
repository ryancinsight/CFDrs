//! Performance metrics, scaling test results, and grading types.

use std::time::Duration;

/// Comprehensive performance metrics for MPI scaling analysis
#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    /// Wall clock time for the entire simulation
    pub total_time: Duration,
    /// Time spent in computation (excluding communication)
    pub computation_time: Duration,
    /// Time spent in MPI communication
    pub communication_time: Duration,
    /// Time spent in load balancing operations
    pub load_balancing_time: Duration,
    /// Time spent in I/O operations
    pub io_time: Duration,
    /// Number of MPI processes used
    pub num_processes: usize,
    /// Problem size (e.g., number of grid cells)
    pub problem_size: usize,
    /// Load imbalance ratio (max load / average load)
    pub load_imbalance_ratio: f64,
    /// Parallel efficiency (ideal_time / actual_time)
    pub parallel_efficiency: f64,
    /// Communication overhead percentage
    pub communication_overhead_percent: f64,
    /// Memory usage per process (MB)
    pub memory_usage_mb: f64,
    /// Network bandwidth utilization (MB/s)
    pub network_bandwidth_mbs: f64,
    /// Cache miss rate (if available)
    pub cache_miss_rate: Option<f64>,
}

impl PerformanceMetrics {
    /// Create new performance metrics with default values
    pub fn new() -> Self {
        Self {
            total_time: Duration::default(),
            computation_time: Duration::default(),
            communication_time: Duration::default(),
            load_balancing_time: Duration::default(),
            io_time: Duration::default(),
            num_processes: 0,
            problem_size: 0,
            load_imbalance_ratio: 1.0,
            parallel_efficiency: 1.0,
            communication_overhead_percent: 0.0,
            memory_usage_mb: 0.0,
            network_bandwidth_mbs: 0.0,
            cache_miss_rate: None,
        }
    }

    /// Calculate derived metrics from raw timing data
    pub fn calculate_derived_metrics(&mut self) {
        if self.total_time.as_secs_f64() > 0.0 {
            self.communication_overhead_percent =
                self.communication_time.as_secs_f64() / self.total_time.as_secs_f64() * 100.0;
        }

        // Parallel efficiency assumes linear scaling from single process
        if self.num_processes > 0 {
            self.parallel_efficiency = 1.0 / self.num_processes as f64;
            if self.total_time.as_secs_f64() > 0.0 {
                // This would need baseline single-process time for accurate calculation
                // For now, use theoretical linear scaling
            }
        }
    }
}

/// Scaling performance grade
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ScalingGrade {
    /// Excellent (≥80% efficiency, ≤10% comm overhead)
    A,
    /// Good (≥70% efficiency, ≤15% comm overhead)
    B,
    /// Acceptable (≥60% efficiency)
    C,
    /// Poor (≥50% efficiency)
    D,
    /// Failing (<50% efficiency)
    F,
}

/// Type of scaling test performed
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalingTestType {
    /// Strong scaling: fixed problem size, increasing cores
    Strong,
    /// Weak scaling: problem size scales with cores
    Weak,
}

/// Results from a scaling benchmark test
#[derive(Debug, Clone)]
pub struct ScalingTestResult {
    /// Core counts tested
    pub core_counts: Vec<usize>,
    /// Performance metrics for each core count
    pub metrics: Vec<PerformanceMetrics>,
    /// Scaling efficiency for each core count (relative to single core)
    pub scaling_efficiency: Vec<f64>,
    /// Communication overhead trend
    pub communication_trend: Vec<f64>,
    /// Load imbalance trend
    pub load_imbalance_trend: Vec<f64>,
    /// Test type (strong or weak scaling)
    pub test_type: ScalingTestType,
    /// Overall assessment
    pub assessment: super::assessment::ScalingAssessment,
}

impl ScalingTestResult {
    /// Calculate scaling efficiency metrics
    pub fn calculate_efficiency(&mut self) {
        if self.core_counts.is_empty() || self.metrics.is_empty() {
            return;
        }

        // Use the first measurement as baseline
        let baseline_time = self.metrics[0].total_time.as_secs_f64();
        let baseline_cores = self.core_counts[0] as f64;

        self.scaling_efficiency = self
            .metrics
            .iter()
            .enumerate()
            .map(|(i, metric)| {
                let cores = self.core_counts[i] as f64;
                let time = metric.total_time.as_secs_f64();

                if time > 0.0 {
                    (baseline_time * baseline_cores) / (time * cores)
                } else {
                    0.0
                }
            })
            .collect();

        self.communication_trend = self
            .metrics
            .iter()
            .map(|m| m.communication_overhead_percent)
            .collect();

        self.load_imbalance_trend = self
            .metrics
            .iter()
            .map(|m| m.load_imbalance_ratio)
            .collect();
    }

    /// Assess overall scaling performance
    pub fn assess_scaling(&mut self) {
        self.assessment =
            super::assessment::ScalingAssessment::assess_from_results(self);
    }
}
