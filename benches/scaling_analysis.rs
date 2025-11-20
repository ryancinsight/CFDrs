//! Scaling analysis framework for CFD operations
//!
//! This module provides comprehensive scaling analysis including:
//! - Weak scaling benchmarks (increasing problem size with cores)
//! - Strong scaling benchmarks (fixed problem size, increasing cores)
//! - Parallel efficiency metrics and analysis
//! - Scaling visualization and reporting

use super::{BenchmarkConfig, PerformanceMetrics};
use criterion::{black_box, BenchmarkId, Criterion, Throughput};
use rayon;
use std::collections::HashMap;
use std::time::{Duration, Instant};

/// Scaling analysis configuration
#[derive(Debug, Clone)]
pub struct ScalingConfig {
    /// Base problem sizes for scaling analysis
    pub base_sizes: Vec<usize>,
    /// Number of threads/processors to test
    pub thread_counts: Vec<usize>,
    /// Scaling type (weak or strong)
    pub scaling_type: ScalingType,
    /// Number of iterations per scaling point
    pub iterations: usize,
    /// Warm-up iterations
    pub warmup_iterations: usize,
}

#[derive(Debug, Clone, Copy)]
pub enum ScalingType {
    /// Weak scaling: problem size increases with thread count
    Weak,
    /// Strong scaling: problem size stays constant
    Strong,
}

impl Default for ScalingConfig {
    fn default() -> Self {
        Self {
            base_sizes: vec![64, 128, 256, 512],
            thread_counts: vec![1, 2, 4, 8, 16],
            scaling_type: ScalingType::Strong,
            iterations: 50,
            warmup_iterations: 5,
        }
    }
}

/// Scaling analysis result
#[derive(Debug, Clone)]
pub struct ScalingResult {
    pub problem_size: usize,
    pub thread_count: usize,
    pub scaling_type: ScalingType,
    pub timing_stats: TimingStatistics,
    pub speedup: f64,
    pub efficiency: f64,
    pub serial_time: Duration,
    pub parallel_time: Duration,
}

/// Timing statistics for scaling analysis
#[derive(Debug, Clone)]
pub struct TimingStatistics {
    pub mean: Duration,
    pub std_dev: Duration,
    pub min: Duration,
    pub max: Duration,
    pub samples: usize,
}

/// Scaling analysis metrics
#[derive(Debug, Clone)]
pub struct ScalingMetrics {
    pub results: Vec<ScalingResult>,
    pub average_speedup: f64,
    pub average_efficiency: f64,
    pub scaling_efficiency: f64,
    pub max_speedup: f64,
    pub optimal_thread_count: usize,
}

/// Run comprehensive scaling analysis
pub fn run_scaling_analysis(config: &ScalingConfig) -> ScalingMetrics {
    println!(
        "Running {} scaling analysis...",
        match config.scaling_type {
            ScalingType::Weak => "weak",
            ScalingType::Strong => "strong",
        }
    );

    let mut results = Vec::new();

    for &base_size in &config.base_sizes {
        println!("Analyzing scaling for base size {}", base_size);

        // Run serial baseline (1 thread)
        let serial_result = benchmark_with_threads(1, base_size, config, config.scaling_type);
        let serial_time = serial_result.timing_stats.mean;

        for &thread_count in &config.thread_counts {
            if thread_count == 1 {
                // Serial result
                results.push(ScalingResult {
                    problem_size: base_size,
                    thread_count: 1,
                    scaling_type: config.scaling_type,
                    timing_stats: serial_result.timing_stats,
                    speedup: 1.0,
                    efficiency: 1.0,
                    serial_time,
                    parallel_time: serial_time,
                });
                continue;
            }

            // Calculate problem size for this thread count
            let problem_size = match config.scaling_type {
                ScalingType::Weak => base_size * thread_count,
                ScalingType::Strong => base_size,
            };

            let parallel_result =
                benchmark_with_threads(thread_count, problem_size, config, config.scaling_type);
            let parallel_time = parallel_result.timing_stats.mean;

            let speedup = serial_time.as_secs_f64() / parallel_time.as_secs_f64();
            let efficiency = speedup / thread_count as f64;

            results.push(ScalingResult {
                problem_size,
                thread_count,
                scaling_type: config.scaling_type,
                timing_stats: parallel_result.timing_stats,
                speedup,
                efficiency,
                serial_time,
                parallel_time,
            });

            println!(
                "  Threads: {}, Size: {}, Speedup: {:.2}x, Efficiency: {:.1}%",
                thread_count,
                problem_size,
                speedup,
                efficiency * 100.0
            );
        }
    }

    // Calculate overall metrics
    let average_speedup = results.iter().map(|r| r.speedup).sum::<f64>() / results.len() as f64;
    let average_efficiency =
        results.iter().map(|r| r.efficiency).sum::<f64>() / results.len() as f64;

    // Calculate scaling efficiency (how well we maintain efficiency as we scale)
    let scaling_efficiency = if results.is_empty() {
        0.0
    } else {
        let max_efficiency = results.iter().map(|r| r.efficiency).fold(0.0, f64::max);
        let min_efficiency = results
            .iter()
            .map(|r| r.efficiency)
            .fold(f64::INFINITY, f64::min);
        if max_efficiency > 0.0 {
            min_efficiency / max_efficiency
        } else {
            0.0
        }
    };

    let max_speedup = results.iter().map(|r| r.speedup).fold(0.0, f64::max);

    let optimal_thread_count = results
        .iter()
        .max_by(|a, b| {
            a.efficiency
                .partial_cmp(&b.efficiency)
                .unwrap_or(std::cmp::Ordering::Equal)
        })
        .map(|r| r.thread_count)
        .unwrap_or(1);

    println!("\n=== Scaling Analysis Summary ===");
    println!("Average speedup: {:.2}x", average_speedup);
    println!("Average efficiency: {:.1}%", average_efficiency * 100.0);
    println!("Scaling efficiency: {:.1}%", scaling_efficiency * 100.0);
    println!("Maximum speedup: {:.2}x", max_speedup);
    println!("Optimal thread count: {}", optimal_thread_count);

    ScalingMetrics {
        results,
        average_speedup,
        average_efficiency,
        scaling_efficiency,
        max_speedup,
        optimal_thread_count,
    }
}

/// Benchmark a CFD operation with specific thread count
fn benchmark_with_threads(
    thread_count: usize,
    problem_size: usize,
    config: &ScalingConfig,
    scaling_type: ScalingType,
) -> BenchmarkResult {
    // Set the number of threads for this benchmark
    let _guard = if thread_count > 1 {
        Some(
            rayon::ThreadPoolBuilder::new()
                .num_threads(thread_count)
                .build_global()
                .expect("Failed to set thread pool"),
        )
    } else {
        None
    };

    // Warm-up phase
    for _ in 0..config.warmup_iterations {
        run_cfd_benchmark_operation(problem_size);
    }

    // Measurement phase
    let mut measurements = Vec::new();

    for _ in 0..config.iterations {
        let start = Instant::now();
        black_box(run_cfd_benchmark_operation(problem_size));
        let elapsed = start.elapsed();
        measurements.push(elapsed);
    }

    // Calculate statistics
    let mean = measurements.iter().sum::<Duration>() / measurements.len() as u32;
    let variance = measurements
        .iter()
        .map(|&d| {
            let diff = if d > mean { d - mean } else { mean - d };
            diff.as_secs_f64().powi(2)
        })
        .sum::<f64>()
        / measurements.len() as f64;

    let std_dev = Duration::from_secs_f64(variance.sqrt());
    let min = *measurements.iter().min().unwrap();
    let max = *measurements.iter().max().unwrap();

    BenchmarkResult {
        timing_stats: TimingStatistics {
            mean,
            std_dev,
            min,
            max,
            samples: measurements.len(),
        },
    }
}

/// Representative CFD benchmark operation
fn run_cfd_benchmark_operation(size: usize) -> f64 {
    use nalgebra::DMatrix;

    // Create CFD simulation fields
    let mut velocity_u = DMatrix::<f64>::zeros(size, size);
    let mut velocity_v = DMatrix::<f64>::zeros(size, size);
    let mut pressure = DMatrix::<f64>::zeros(size, size);

    // Simulate CFD computation (pressure Poisson equation solver)
    let mut residual = 0.0;
    let iterations = size.min(20); // Limit iterations for scaling analysis

    for _ in 0..iterations {
        let mut max_residual = 0.0;

        // Parallel computation using rayon
        use rayon::prelude::*;
        velocity_u.par_iter_mut().enumerate().for_each(|(i, row)| {
            row.iter_mut().enumerate().for_each(|(j, cell)| {
                if i > 0 && i < size - 1 && j > 0 && j < size - 1 {
                    // Simple Jacobi iteration for pressure
                    let new_value = (pressure[(i - 1, j)]
                        + pressure[(i + 1, j)]
                        + pressure[(i, j - 1)]
                        + pressure[(i, j + 1)])
                        / 4.0;
                    let residual = (new_value - *cell).abs();
                    *cell = new_value;
                    // Note: This is a simplified version - in practice we'd use atomic operations
                }
            });
        });

        residual = max_residual;
    }

    residual
}

/// Result of a single benchmark run
#[derive(Debug)]
struct BenchmarkResult {
    timing_stats: TimingStatistics,
}

/// Benchmark scaling behavior with criterion
pub fn benchmark_scaling_behavior(c: &mut Criterion, config: &BenchmarkConfig) {
    let scaling_config = ScalingConfig::default();

    let mut group = c.benchmark_group("scaling_analysis");

    for &base_size in &scaling_config.base_sizes {
        for &thread_count in &scaling_config.thread_counts {
            let problem_size = match scaling_config.scaling_type {
                ScalingType::Weak => base_size * thread_count,
                ScalingType::Strong => base_size,
            };

            group.throughput(Throughput::Elements((problem_size * problem_size) as u64));

            group.bench_with_input(
                BenchmarkId::new(
                    format!(
                        "scaling_{}_{}threads",
                        match scaling_config.scaling_type {
                            ScalingType::Weak => "weak",
                            ScalingType::Strong => "strong",
                        },
                        thread_count
                    ),
                    problem_size,
                ),
                &(thread_count, problem_size),
                |b, &(thread_count, problem_size)| {
                    // Set thread count for this benchmark
                    let _guard = if thread_count > 1 {
                        Some(
                            rayon::ThreadPoolBuilder::new()
                                .num_threads(thread_count)
                                .build_global()
                                .expect("Failed to set thread pool"),
                        )
                    } else {
                        None
                    };

                    b.iter(|| {
                        black_box(run_cfd_benchmark_operation(problem_size));
                    });
                },
            );
        }
    }

    group.finish();

    // Run comprehensive scaling analysis
    let metrics = run_scaling_analysis(&scaling_config);
    save_scaling_results(&metrics);
}

/// Save scaling analysis results
fn save_scaling_results(metrics: &ScalingMetrics) {
    use std::fs;
    use std::path::Path;

    let output_dir = "performance_results/scaling";
    let dir_path = Path::new(output_dir);

    if let Err(e) = fs::create_dir_all(dir_path) {
        eprintln!(
            "Failed to create scaling output directory {}: {}",
            output_dir, e
        );
        return;
    }

    // Save detailed results
    let results_path = dir_path.join("scaling_results.json");
    if let Ok(json) = serde_json::to_string_pretty(&metrics.results) {
        if let Err(e) = fs::write(&results_path, json) {
            eprintln!(
                "Failed to save scaling results to {}: {}",
                results_path.display(),
                e
            );
        } else {
            println!("Scaling results saved to {}", results_path.display());
        }
    }

    // Save summary metrics
    let summary = serde_json::json!({
        "average_speedup": metrics.average_speedup,
        "average_efficiency": metrics.average_efficiency,
        "scaling_efficiency": metrics.scaling_efficiency,
        "max_speedup": metrics.max_speedup,
        "optimal_thread_count": metrics.optimal_thread_count,
    });

    let summary_path = dir_path.join("scaling_summary.json");
    if let Ok(json) = serde_json::to_string_pretty(&summary) {
        if let Err(e) = fs::write(&summary_path, json) {
            eprintln!(
                "Failed to save scaling summary to {}: {}",
                summary_path.display(),
                e
            );
        } else {
            println!("Scaling summary saved to {}", summary_path.display());
        }
    }
}

/// Generate scaling recommendations based on analysis
pub fn generate_scaling_recommendations(metrics: &ScalingMetrics) -> Vec<String> {
    let mut recommendations = Vec::new();

    if metrics.average_efficiency < 0.5 {
        recommendations.push(format!(
            "Poor parallel efficiency ({:.1}%). Consider optimizing parallel algorithms or reducing communication overhead.",
            metrics.average_efficiency * 100.0
        ));
    }

    if metrics.scaling_efficiency < 0.7 {
        recommendations.push(format!(
            "Poor scaling efficiency ({:.1}%). Performance degrades significantly with increased thread count.",
            metrics.scaling_efficiency * 100.0
        ));
    }

    if metrics.optimal_thread_count > 1 {
        recommendations.push(format!(
            "Optimal thread count is {}. Using more threads may decrease efficiency.",
            metrics.optimal_thread_count
        ));
    } else {
        recommendations.push(
            "Single-threaded performance is optimal. Consider if parallelization is necessary."
                .to_string(),
        );
    }

    if metrics.max_speedup > 10.0 {
        recommendations.push(format!(
            "Excellent speedup achieved ({:.1}x). Parallel implementation is highly effective.",
            metrics.max_speedup
        ));
    } else if metrics.max_speedup < 2.0 {
        recommendations.push(format!(
            "Limited speedup ({:.1}x). Consider algorithmic improvements or different parallelization strategy.",
            metrics.max_speedup
        ));
    }

    recommendations
}
