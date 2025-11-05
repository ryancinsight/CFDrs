//! Comprehensive CFD Performance Benchmarking Suite
//!
//! This benchmark suite provides comprehensive performance analysis for all major CFD operations,
//! including memory profiling, scaling analysis, and regression detection capabilities.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use std::time::{Duration, Instant};
use std::collections::HashMap;

// Core CFD operation benchmarks
mod core_operations;
mod memory_profiling;
mod scaling_analysis;
mod regression_detection;

// Re-export for external use
pub use core_operations::*;
pub use memory_profiling::*;
pub use scaling_analysis::*;
pub use regression_detection::*;

/// Main benchmark configuration
#[derive(Debug, Clone)]
pub struct BenchmarkConfig {
    pub problem_sizes: Vec<usize>,
    pub iterations: usize,
    pub warmup_iterations: usize,
    pub enable_memory_profiling: bool,
    pub enable_scaling_analysis: bool,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            problem_sizes: vec![32, 64, 128, 256, 512],
            iterations: 100,
            warmup_iterations: 10,
            enable_memory_profiling: true,
            enable_scaling_analysis: true,
        }
    }
}

/// Performance metrics collector
#[derive(Debug, Clone)]
pub struct PerformanceMetrics {
    pub operation_name: String,
    pub problem_size: usize,
    pub execution_time: Duration,
    pub memory_usage: Option<u64>,
    pub throughput: Option<f64>,
    pub scaling_efficiency: Option<f64>,
}

impl PerformanceMetrics {
    pub fn new(operation_name: String, problem_size: usize) -> Self {
        Self {
            operation_name,
            problem_size,
            execution_time: Duration::default(),
            memory_usage: None,
            throughput: None,
            scaling_efficiency: None,
        }
    }

    pub fn with_timing(mut self, time: Duration) -> Self {
        self.execution_time = time;
        self
    }

    pub fn with_memory(mut self, memory: u64) -> Self {
        self.memory_usage = Some(memory);
        self
    }

    pub fn with_throughput(mut self, throughput: f64) -> Self {
        self.throughput = Some(throughput);
        self
    }

    pub fn with_scaling_efficiency(mut self, efficiency: f64) -> Self {
        self.scaling_efficiency = Some(efficiency);
        self
    }
}

/// Benchmark runner with comprehensive analysis
pub struct ComprehensiveBenchmarkRunner {
    config: BenchmarkConfig,
    results: Vec<PerformanceMetrics>,
}

impl ComprehensiveBenchmarkRunner {
    pub fn new(config: BenchmarkConfig) -> Self {
        Self {
            config,
            results: Vec::new(),
        }
    }

    pub fn run_comprehensive_suite(&mut self, c: &mut Criterion) {
        // Core CFD operations
        self.benchmark_grid_operations(c);
        self.benchmark_solver_operations(c);
        self.benchmark_turbulence_operations(c);
        self.benchmark_gpu_operations(c);

        // Advanced analysis
        if self.config.enable_memory_profiling {
            self.benchmark_memory_usage(c);
        }

        if self.config.enable_scaling_analysis {
            self.benchmark_parallel_scaling(c);
        }

        // Regression detection
        self.perform_regression_analysis();
    }

    fn benchmark_grid_operations(&mut self, c: &mut Criterion) {
        core_operations::benchmark_grid_operations(c, &self.config);
    }

    fn benchmark_solver_operations(&mut self, c: &mut Criterion) {
        core_operations::benchmark_solver_operations(c, &self.config);
    }

    fn benchmark_turbulence_operations(&mut self, c: &mut Criterion) {
        core_operations::benchmark_turbulence_operations(c, &self.config);
    }

    fn benchmark_gpu_operations(&mut self, c: &mut Criterion) {
        core_operations::benchmark_gpu_operations(c, &self.config);
    }

    fn benchmark_memory_usage(&mut self, c: &mut Criterion) {
        memory_profiling::benchmark_memory_usage(c, &self.config);
    }

    fn benchmark_parallel_scaling(&mut self, c: &mut Criterion) {
        scaling_analysis::benchmark_parallel_scaling(c, &self.config);
    }

    fn perform_regression_analysis(&mut self) {
        regression_detection::analyze_performance_regression(&self.results);
    }

    pub fn get_results(&self) -> &[PerformanceMetrics] {
        &self.results
    }
}

// Main benchmark functions
fn benchmark_comprehensive_cfd_operations(c: &mut Criterion) {
    let config = BenchmarkConfig::default();
    let mut runner = ComprehensiveBenchmarkRunner::new(config);
    runner.run_comprehensive_suite(c);
}

fn benchmark_memory_profiling_suite(c: &mut Criterion) {
    let config = BenchmarkConfig {
        enable_memory_profiling: true,
        enable_scaling_analysis: false,
        ..Default::default()
    };
    let mut runner = ComprehensiveBenchmarkRunner::new(config);
    runner.benchmark_memory_usage(c);
}

fn benchmark_scaling_analysis_suite(c: &mut Criterion) {
    let config = BenchmarkConfig {
        enable_memory_profiling: false,
        enable_scaling_analysis: true,
        ..Default::default()
    };
    let mut runner = ComprehensiveBenchmarkRunner::new(config);
    runner.benchmark_parallel_scaling(c);
}

criterion_group!(
    benches,
    benchmark_comprehensive_cfd_operations,
    benchmark_memory_profiling_suite,
    benchmark_scaling_analysis_suite
);

criterion_main!(benches);
