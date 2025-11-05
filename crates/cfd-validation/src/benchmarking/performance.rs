//! Performance benchmarking for CFD operations
//!
//! Measures execution time, throughput, and identifies performance bottlenecks
//! in CFD algorithms and data structures.

use super::utils::{BenchmarkStats, BenchmarkTimer};
use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Scalar};
use nalgebra_sparse::CsrMatrix;
use std::collections::HashMap;
use std::fmt;

/// Result of a single performance benchmark
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TimingResult {
    /// Name of the operation being benchmarked
    pub operation_name: String,
    /// Individual timing measurements (in seconds)
    pub measurements: Vec<f64>,
    /// Statistical summary
    pub stats: BenchmarkStats,
    /// Additional metadata
    pub metadata: HashMap<String, String>,
}

impl TimingResult {
    /// Create a new timing result
    pub fn new(operation_name: String, measurements: Vec<f64>) -> Self {
        let stats = BenchmarkStats::from_measurements(&measurements);
        Self {
            operation_name,
            measurements,
            stats,
            metadata: HashMap::new(),
        }
    }

    /// Add metadata
    pub fn with_metadata(mut self, key: String, value: String) -> Self {
        self.metadata.insert(key, value);
        self
    }

    /// Check if benchmark is stable (low coefficient of variation)
    pub fn is_stable(&self, max_cv: f64) -> bool {
        self.stats.is_stable(max_cv)
    }
}

impl fmt::Display for TimingResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {:.3}±{:.3}ms ({} samples)",
               self.operation_name,
               self.stats.mean * 1000.0,
               self.stats.std_dev * 1000.0,
               self.stats.samples)?;

        if !self.metadata.is_empty() {
            write!(f, " [")?;
            for (i, (key, value)) in self.metadata.iter().enumerate() {
                if i > 0 { write!(f, ", ")?; }
                write!(f, "{}={}", key, value)?;
            }
            write!(f, "]")?;
        }

        Ok(())
    }
}

/// Performance benchmark runner
pub struct PerformanceBenchmark {
    /// Number of warm-up iterations
    warm_up_iterations: usize,
    /// Number of measurement iterations
    measurement_iterations: usize,
    /// Minimum measurement time per iteration (seconds)
    min_measurement_time: f64,
    /// Maximum coefficient of variation for stable results
    max_cv_threshold: f64,
}

impl PerformanceBenchmark {
    /// Create a new performance benchmark with default settings
    pub fn new() -> Self {
        Self {
            warm_up_iterations: 5,
            measurement_iterations: 10,
            min_measurement_time: 0.001, // 1ms minimum
            max_cv_threshold: 0.05, // 5% CV threshold
        }
    }

    /// Configure warm-up iterations
    pub fn with_warm_up(mut self, iterations: usize) -> Self {
        self.warm_up_iterations = iterations;
        self
    }

    /// Configure measurement iterations
    pub fn with_measurements(mut self, iterations: usize) -> Self {
        self.measurement_iterations = iterations;
        self
    }

    /// Configure minimum measurement time
    pub fn with_min_time(mut self, seconds: f64) -> Self {
        self.min_measurement_time = seconds;
        self
    }

    /// Configure stability threshold
    pub fn with_stability_threshold(mut self, cv: f64) -> Self {
        self.max_cv_threshold = cv;
        self
    }

    /// Benchmark a closure that returns a result
    pub fn benchmark<F, T, E>(&self, operation_name: &str, mut operation: F) -> Result<TimingResult>
    where
        F: FnMut() -> std::result::Result<T, E>,
        E: std::fmt::Display,
    {
        // Warm-up phase
        for _ in 0..self.warm_up_iterations {
            operation().map_err(|e| Error::InvalidInput(format!("Warm-up failed: {}", e)))?;
        }

        // Measurement phase
        let mut measurements = Vec::with_capacity(self.measurement_iterations);
        let mut timer = BenchmarkTimer::new();

        for _ in 0..self.measurement_iterations {
            timer.start();

            // Ensure minimum measurement time
            let start_instant = std::time::Instant::now();
            operation().map_err(|e| Error::InvalidInput(format!("Operation failed: {}", e)))?;

            // Pad measurement if too fast
            let elapsed = start_instant.elapsed().as_secs_f64();
            if elapsed < self.min_measurement_time {
                let remaining = self.min_measurement_time - elapsed;
                std::thread::sleep(std::time::Duration::from_secs_f64(remaining));
            }

            let duration = timer.stop();
            measurements.push(duration.as_secs_f64());
        }

        Ok(TimingResult::new(operation_name.to_string(), measurements))
    }

    /// Benchmark a simple closure (no result)
    pub fn benchmark_simple<F>(&self, operation_name: &str, mut operation: F) -> Result<TimingResult>
    where
        F: FnMut(),
    {
        self.benchmark(operation_name, || {
            operation();
            Ok::<(), std::convert::Infallible>(())
        })
    }

    /// Run multiple benchmarks and return results
    pub fn benchmark_suite(&self, benchmarks: Vec<(&str, Box<dyn FnMut()>)>) -> Result<Vec<TimingResult>> {
        let mut results = Vec::new();

        for (name, mut operation) in benchmarks {
            let result = self.benchmark_simple(name, &mut operation)?;
            results.push(result);
        }

        Ok(results)
    }
}

impl Default for PerformanceBenchmark {
    fn default() -> Self {
        Self::new()
    }
}

/// CFD-specific performance benchmarks
pub struct CfdPerformanceBenchmarks {
    benchmark: PerformanceBenchmark,
}

impl CfdPerformanceBenchmarks {
    pub fn new() -> Self {
        Self {
            benchmark: PerformanceBenchmark::new()
                .with_warm_up(3)
                .with_measurements(20)
                .with_min_time(0.01) // 10ms minimum for CFD operations
                .with_stability_threshold(0.10), // 10% CV for CFD (more variable)
        }
    }

    /// Benchmark matrix-vector multiplication (key CFD operation)
    pub fn benchmark_spmv<T>(&self, matrix: &CsrMatrix<T>, vector: &[T]) -> Result<TimingResult>
    where
        T: RealField + Copy + Scalar + std::fmt::Display,
    {
        let mut result = vec![T::zero(); matrix.nrows()];

        self.benchmark.benchmark_simple(
            &format!("SPMV_{}x{}", matrix.nrows(), matrix.ncols()),
            || {
                let matrix = &matrix;
                let vector_dv = nalgebra::DVector::from_vec(vector.to_vec());
                let mut result_dv = nalgebra::DVector::from_vec(result.to_vec());
                cfd_math::sparse::spmv(matrix, &vector_dv, &mut result_dv);
            }
        )
    }

    /// Benchmark linear solver convergence
    pub fn benchmark_solver_convergence<T>(
        &self,
        solver_name: &str,
        solver_fn: impl Fn() -> Result<Vec<T>>
    ) -> Result<TimingResult>
    where
        T: nalgebra::RealField + Copy,
    {
        self.benchmark.benchmark(
            &format!("Solver_{}", solver_name),
            || solver_fn()
        )
    }

    /// Benchmark grid operations
    pub fn benchmark_grid_operation<T>(
        &self,
        operation_name: &str,
        operation: impl Fn()
    ) -> Result<TimingResult> {
        self.benchmark.benchmark_simple(
            &format!("Grid_{}", operation_name),
            operation
        )
    }

    /// Benchmark CFD time stepping
    pub fn benchmark_time_step<T>(
        &self,
        scheme_name: &str,
        time_step_fn: impl Fn() -> Result<T>
    ) -> Result<TimingResult> {
        self.benchmark.benchmark(
            &format!("TimeStep_{}", scheme_name),
            time_step_fn
        )
    }

    /// Run comprehensive CFD benchmark suite
    pub fn run_cfd_suite(&self) -> Result<Vec<TimingResult>> {
        println!("Running CFD Performance Benchmark Suite...");

        let mut benchmarks = Vec::new();

        // Matrix operations benchmark
        benchmarks.push((
            "Matrix_Assembly_64x64",
            Box::new(|| {
                let mut builder = cfd_math::sparse::SparseMatrixBuilder::new(64 * 64, 64 * 64);
                for i in 0..64 {
                    for j in 0..64 {
                        let idx = i * 64 + j;
                        if i > 0 { builder.add_entry(idx, idx - 64, -1.0).unwrap(); }
                        if i < 63 { builder.add_entry(idx, idx + 64, -1.0).unwrap(); }
                        if j > 0 { builder.add_entry(idx, idx - 1, -1.0).unwrap(); }
                        if j < 63 { builder.add_entry(idx, idx + 1, -1.0).unwrap(); }
                        builder.add_entry(idx, idx, 4.0).unwrap();
                    }
                }
                let _matrix = builder.build().unwrap();
            }) as Box<dyn FnMut()>
        ));

        // Vector operations benchmark
        benchmarks.push((
            "Vector_Operations_1000",
            Box::new(|| {
                let mut v1 = vec![1.0f64; 1000];
                let v2 = vec![2.0f64; 1000];
                for i in 0..1000 {
                    v1[i] = v1[i] * v2[i] + v1[i].sin();
                }
            }) as Box<dyn FnMut()>
        ));

        // Memory allocation benchmark
        benchmarks.push((
            "Memory_Allocation_1M",
            Box::new(|| {
                let _data = vec![0.0f64; 1_000_000];
            }) as Box<dyn FnMut()>
        ));

        self.benchmark.benchmark_suite(benchmarks)
    }
}

impl Default for CfdPerformanceBenchmarks {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_performance_benchmark() {
        let benchmark = PerformanceBenchmark::new();

        let result = benchmark.benchmark_simple("test_operation", || {
            std::thread::sleep(std::time::Duration::from_millis(1));
        }).unwrap();

        assert_eq!(result.operation_name, "test_operation");
        assert!(!result.measurements.is_empty());
        assert!(result.stats.mean > 0.0);
        assert!(result.is_stable(0.5)); // Should be reasonably stable
    }

    #[test]
    fn test_cfd_benchmarks() {
        let cfd_benchmarks = CfdPerformanceBenchmarks::new();
        let results = cfd_benchmarks.run_cfd_suite().unwrap();

        assert!(!results.is_empty());
        for result in results {
            assert!(result.stats.mean > 0.0);
            assert!(result.stats.samples > 0);
            println!("{}", result);
        }
    }

    #[test]
    fn test_timing_result_display() {
        let measurements = vec![0.1, 0.105, 0.095];
        let result = TimingResult::new("test".to_string(), measurements)
            .with_metadata("grid".to_string(), "64x64".to_string());

        let display = format!("{}", result);
        assert!(display.contains("test"));
        assert!(display.contains("grid=64x64"));
    }
}
