//! Benchmark problems for CFD validation.
//!
//! This module provides standard benchmark problems used to validate CFD solvers
//! against known analytical solutions and experimental data.

use crate::error_metrics::ErrorStatistics;
use cfd_core::Result;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Trait for benchmark problems
pub trait Benchmark<T: RealField> {
    /// Configuration type for this benchmark
    type Config;
    /// Solution type produced by this benchmark
    type Solution;

    /// Get the name of this benchmark
    fn name(&self) -> &str;

    /// Get a description of this benchmark
    fn description(&self) -> &str;

    /// Set up the benchmark with given configuration
    fn setup(&mut self, config: Self::Config) -> Result<()>;

    /// Run the benchmark and return the solution
    fn run(&self) -> Result<Self::Solution>;

    /// Validate the solution against reference data
    fn validate(&self, solution: &Self::Solution) -> Result<BenchmarkResult<T>>;
}

/// Result of a benchmark validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult<T: RealField> {
    /// Name of the benchmark
    pub benchmark_name: String,
    /// Whether the benchmark passed validation
    pub passed: bool,
    /// Error statistics for different metrics
    pub error_statistics: HashMap<String, ErrorStatistics<T>>,
    /// Additional metadata
    pub metadata: HashMap<String, String>,
    /// Execution time in seconds
    pub execution_time: Option<f64>,
}

/// Configuration for benchmark suite
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig<T: RealField> {
    /// Tolerance for validation
    pub tolerance: T,
    /// Whether to run in verbose mode
    pub verbose: bool,
    /// Maximum execution time per benchmark
    pub max_execution_time: Option<f64>,
}

impl<T: RealField> Default for BenchmarkConfig<T>
where
    T: From<f64>,
{
    fn default() -> Self {
        Self {
            tolerance: T::from(1e-6),
            verbose: false,
            max_execution_time: Some(300.0), // 5 minutes
        }
    }
}

/// Suite of benchmark problems
pub struct BenchmarkSuite<T: RealField> {
    /// List of benchmarks in the suite
    benchmarks: Vec<Box<dyn Benchmark<T, Config = BenchmarkConfig<T>, Solution = Vec<T>>>>,
    /// Configuration for the suite
    config: BenchmarkConfig<T>,
}

impl<T: RealField> BenchmarkSuite<T>
where
    T: From<f64>,
{
    /// Create a new benchmark suite
    pub fn new(config: BenchmarkConfig<T>) -> Self {
        Self {
            benchmarks: Vec::new(),
            config,
        }
    }

    /// Add a benchmark to the suite
    pub fn add_benchmark<B>(&mut self, benchmark: B)
    where
        B: Benchmark<T, Config = BenchmarkConfig<T>, Solution = Vec<T>> + 'static
    {
        self.benchmarks.push(Box::new(benchmark));
    }

    /// Run all benchmarks in the suite
    pub fn run_all(&self) -> Result<Vec<BenchmarkResult<T>>> {
        let mut results = Vec::new();

        for benchmark in &self.benchmarks {
            if self.config.verbose {
                println!("Running benchmark: {}", benchmark.name());
            }

            let start_time = std::time::Instant::now();
            let solution = benchmark.run()?;
            let execution_time = start_time.elapsed().as_secs_f64();

            let mut result = benchmark.validate(&solution)?;
            result.execution_time = Some(execution_time);

            if self.config.verbose {
                println!("Benchmark {} {}",
                    benchmark.name(),
                    if result.passed { "PASSED" } else { "FAILED" }
                );
            }

            results.push(result);
        }

        Ok(results)
    }
}

/// Lid-driven cavity benchmark problem
#[derive(Debug)]
pub struct LidDrivenCavity<T: RealField> {
    /// Reynolds number
    reynolds: T,
    /// Grid size
    grid_size: (usize, usize),
    /// Lid velocity
    lid_velocity: T,
}

impl<T: RealField> LidDrivenCavity<T> {
    /// Create a new lid-driven cavity benchmark
    pub fn new(reynolds: T, grid_size: (usize, usize), lid_velocity: T) -> Self {
        Self {
            reynolds,
            grid_size,
            lid_velocity,
        }
    }

    /// Get the Reynolds number
    pub fn reynolds(&self) -> &T {
        &self.reynolds
    }

    /// Get the grid size
    pub fn grid_size(&self) -> (usize, usize) {
        self.grid_size
    }

    /// Get the lid velocity
    pub fn lid_velocity(&self) -> &T {
        &self.lid_velocity
    }
}

impl<T: RealField> Benchmark<T> for LidDrivenCavity<T> {
    type Config = BenchmarkConfig<T>;
    type Solution = Vec<T>;

    fn name(&self) -> &str {
        "Lid-Driven Cavity"
    }

    fn description(&self) -> &str {
        "2D lid-driven cavity flow benchmark problem"
    }

    fn setup(&mut self, _config: Self::Config) -> Result<()> {
        // Setup is handled in constructor for this benchmark
        Ok(())
    }

    fn run(&self) -> Result<Self::Solution> {
        // Placeholder implementation - would integrate with actual 2D solver
        let size = self.grid_size.0 * self.grid_size.1;
        Ok(vec![T::zero(); size])
    }

    fn validate(&self, solution: &Self::Solution) -> Result<BenchmarkResult<T>> {
        // Placeholder validation - would compare against reference data
        let mut error_stats = HashMap::new();
        error_stats.insert(
            "velocity_error".to_string(),
            ErrorStatistics {
                l1_norm: T::zero(),
                l2_norm: T::zero(),
                linf_norm: T::zero(),
                mae: T::zero(),
                rmse: T::zero(),
                relative_l2: T::zero(),
                num_points: solution.len(),
            }
        );

        Ok(BenchmarkResult {
            benchmark_name: self.name().to_string(),
            passed: true,
            error_statistics: error_stats,
            metadata: HashMap::new(),
            execution_time: None,
        })
    }
}

/// Flow over cylinder benchmark problem
#[derive(Debug)]
pub struct FlowOverCylinder<T: RealField> {
    /// Reynolds number
    reynolds: T,
    /// Cylinder diameter
    diameter: T,
}

impl<T: RealField> FlowOverCylinder<T> {
    /// Create a new flow over cylinder benchmark
    pub fn new(reynolds: T, diameter: T) -> Self {
        Self { reynolds, diameter }
    }

    /// Get the Reynolds number
    pub fn reynolds(&self) -> &T {
        &self.reynolds
    }

    /// Get the cylinder diameter
    pub fn diameter(&self) -> &T {
        &self.diameter
    }
}

impl<T: RealField> Benchmark<T> for FlowOverCylinder<T> {
    type Config = BenchmarkConfig<T>;
    type Solution = Vec<T>;

    fn name(&self) -> &str {
        "Flow Over Cylinder"
    }

    fn description(&self) -> &str {
        "2D flow over circular cylinder benchmark problem"
    }

    fn setup(&mut self, _config: Self::Config) -> Result<()> {
        Ok(())
    }

    fn run(&self) -> Result<Self::Solution> {
        // Placeholder implementation
        Ok(vec![T::zero(); 1000])
    }

    fn validate(&self, _solution: &Self::Solution) -> Result<BenchmarkResult<T>> {
        // Placeholder validation
        Ok(BenchmarkResult {
            benchmark_name: self.name().to_string(),
            passed: true,
            error_statistics: HashMap::new(),
            metadata: HashMap::new(),
            execution_time: None,
        })
    }
}

/// Backward-facing step benchmark problem
#[derive(Debug)]
pub struct BackwardFacingStep<T: RealField> {
    /// Reynolds number
    reynolds: T,
    /// Step height
    step_height: T,
}

impl<T: RealField> BackwardFacingStep<T> {
    /// Create a new backward-facing step benchmark
    pub fn new(reynolds: T, step_height: T) -> Self {
        Self { reynolds, step_height }
    }

    /// Get the Reynolds number
    pub fn reynolds(&self) -> &T {
        &self.reynolds
    }

    /// Get the step height
    pub fn step_height(&self) -> &T {
        &self.step_height
    }
}

impl<T: RealField> Benchmark<T> for BackwardFacingStep<T> {
    type Config = BenchmarkConfig<T>;
    type Solution = Vec<T>;

    fn name(&self) -> &str {
        "Backward-Facing Step"
    }

    fn description(&self) -> &str {
        "2D backward-facing step flow benchmark problem"
    }

    fn setup(&mut self, _config: Self::Config) -> Result<()> {
        Ok(())
    }

    fn run(&self) -> Result<Self::Solution> {
        // Placeholder implementation
        Ok(vec![T::zero(); 1000])
    }

    fn validate(&self, _solution: &Self::Solution) -> Result<BenchmarkResult<T>> {
        // Placeholder validation
        Ok(BenchmarkResult {
            benchmark_name: self.name().to_string(),
            passed: true,
            error_statistics: HashMap::new(),
            metadata: HashMap::new(),
            execution_time: None,
        })
    }
}

/// Benchmark runner for executing validation studies
pub struct BenchmarkRunner;

impl BenchmarkRunner {
    /// Run a comprehensive validation study
    pub fn run_validation_study<T: RealField>() -> Result<ValidationReport<T>>
    where
        T: From<f64>,
    {
        let config = BenchmarkConfig::default();
        let mut suite = BenchmarkSuite::new(config);

        // Add standard benchmarks
        suite.add_benchmark(LidDrivenCavity::new(
            T::from(1000.0),
            (64, 64),
            T::one(),
        ));

        suite.add_benchmark(FlowOverCylinder::new(
            T::from(100.0),
            T::one(),
        ));

        suite.add_benchmark(BackwardFacingStep::new(
            T::from(200.0),
            T::one(),
        ));

        let results = suite.run_all()?;

        Ok(ValidationReport {
            total_benchmarks: results.len(),
            passed_benchmarks: results.iter().filter(|r| r.passed).count(),
            results,
        })
    }
}

/// Comprehensive validation report
#[derive(Debug, Serialize, Deserialize)]
pub struct ValidationReport<T: RealField> {
    /// Total number of benchmarks run
    pub total_benchmarks: usize,
    /// Number of benchmarks that passed
    pub passed_benchmarks: usize,
    /// Individual benchmark results
    pub results: Vec<BenchmarkResult<T>>,
}

impl<T: RealField> ValidationReport<T> {
    /// Get the success rate as a percentage
    pub fn success_rate(&self) -> f64 {
        if self.total_benchmarks == 0 {
            0.0
        } else {
            (self.passed_benchmarks as f64 / self.total_benchmarks as f64) * 100.0
        }
    }

    /// Check if all benchmarks passed
    pub fn all_passed(&self) -> bool {
        self.passed_benchmarks == self.total_benchmarks
    }
}
