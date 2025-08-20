//! CFD benchmark suite for validation
//!
//! This module provides standard CFD benchmark problems for validating
//! solver implementations against known solutions.

use nalgebra::RealField;
use cfd_core::Result;
use serde::{Deserialize, Serialize};

pub mod cavity;
pub mod cylinder;
pub mod step;
pub mod runner;

pub use cavity::LidDrivenCavity;
pub use cylinder::FlowOverCylinder;
pub use step::BackwardFacingStep;
pub use runner::{BenchmarkRunner, ValidationReport};

/// Trait for CFD benchmark problems
pub trait Benchmark<T: RealField> {
    /// Get the name of the benchmark
    fn name(&self) -> &str;
    
    /// Get a description of the benchmark problem
    fn description(&self) -> &str;
    
    /// Run the benchmark and return results
    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>>;
    
    /// Get the reference solution (if available)
    fn reference_solution(&self) -> Option<BenchmarkResult<T>>;
    
    /// Validate results against reference solution
    fn validate(&self, result: &BenchmarkResult<T>) -> Result<bool>;
}

/// Results from running a benchmark
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkResult<T: RealField> {
    /// Name of the benchmark
    pub name: String,
    /// Computed solution values
    pub values: Vec<T>,
    /// Error metrics
    pub errors: Vec<T>,
    /// Convergence history
    pub convergence: Vec<T>,
    /// Execution time in seconds
    pub execution_time: f64,
    /// Additional metadata
    pub metadata: std::collections::HashMap<String, String>,
}

/// Configuration for running benchmarks
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BenchmarkConfig<T: RealField> {
    /// Grid resolution
    pub resolution: usize,
    /// Convergence tolerance
    pub tolerance: T,
    /// Maximum iterations
    pub max_iterations: usize,
    /// Reynolds number
    pub reynolds_number: T,
    /// Time step (for transient problems)
    pub time_step: Option<T>,
    /// Enable parallel execution
    pub parallel: bool,
}

impl<T: RealField> Default for BenchmarkConfig<T> {
    fn default() -> Self {
        Self {
            resolution: 64,
            tolerance: T::from_f64(1e-6).unwrap_or_else(|| T::from_f64(0.000001).unwrap()),
            max_iterations: 1000,
            reynolds_number: T::from_f64(100.0).unwrap_or_else(|| T::from_i32(100).unwrap()),
            time_step: None,
            parallel: true,
        }
    }
}

/// Benchmark suite for running multiple benchmarks
pub struct BenchmarkSuite<T: RealField> {
    benchmarks: Vec<Box<dyn Benchmark<T>>>,
}

impl<T: RealField> BenchmarkSuite<T> {
    /// Create a new benchmark suite
    pub fn new() -> Self {
        Self {
            benchmarks: Vec::new(),
        }
    }
    
    /// Add a benchmark to the suite
    pub fn add_benchmark(&mut self, benchmark: Box<dyn Benchmark<T>>) {
        self.benchmarks.push(benchmark);
    }
    
    /// Run all benchmarks
    pub fn run_all(&self, config: &BenchmarkConfig<T>) -> Vec<Result<BenchmarkResult<T>>> {
        self.benchmarks
            .iter()
            .map(|b| b.run(config))
            .collect()
    }
}

impl<T: RealField> Default for BenchmarkSuite<T> {
    fn default() -> Self {
        Self::new()
    }
}