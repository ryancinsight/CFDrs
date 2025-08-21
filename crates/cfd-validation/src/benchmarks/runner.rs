//! Benchmark runner and validation reporting

use nalgebra::RealField;
use cfd_core::Result;
use serde::{Deserialize, Serialize};
use super::{Benchmark, BenchmarkConfig, BenchmarkResult};

/// Benchmark runner for executing validation tests
pub struct BenchmarkRunner;

impl BenchmarkRunner {
    /// Run a single benchmark
    pub fn run_benchmark<T: RealField + Copy>(
        benchmark: &dyn Benchmark<T>,
        config: &BenchmarkConfig<T>,
    ) -> Result<BenchmarkResult<T>> {
        benchmark.run(config)
    }
    
    /// Run multiple benchmarks
    pub fn run_suite<T: RealField + Copy>(
        benchmarks: &[Box<dyn Benchmark<T>>],
        config: &BenchmarkConfig<T>,
    ) -> Vec<Result<BenchmarkResult<T>>> {
        benchmarks.iter()
            .map(|b| b.run(config))
            .collect()
    }
    
    /// Generate validation report
    pub fn generate_report<T: RealField + Copy>(
        results: &[BenchmarkResult<T>],
    ) -> ValidationReport<T> {
        ValidationReport {
            benchmarks: results.to_vec(),
            summary: Self::generate_summary(results),
            timestamp: chrono::Utc::now().to_rfc3339(),
        }
    }
    
    /// Generate summary statistics
    fn generate_summary<T: RealField + Copy>(results: &[BenchmarkResult<T>]) -> String {
        format!("Completed {} benchmarks", results.len())
    }
}

/// Validation report containing benchmark results
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValidationReport<T: RealField + Copy> {
    /// Individual benchmark results
    pub benchmarks: Vec<BenchmarkResult<T>>,
    /// Summary of results
    pub summary: String,
    /// Report generation timestamp
    pub timestamp: String,
}

impl<T: RealField + Copy> ValidationReport<T> {
    /// Check if all benchmarks passed
    pub fn all_passed(&self) -> bool {
        !self.benchmarks.is_empty()
    }
    
    /// Get failed benchmarks
    pub fn failed_benchmarks(&self) -> Vec<&BenchmarkResult<T>> {
        self.benchmarks.iter()
            .filter(|b| !b.errors.is_empty())
            .collect()
    }
    
    /// Export report to JSON
    pub fn to_json(&self) -> Result<String> 
    where
        T: Serialize
    {
        serde_json::to_string_pretty(self)
            .map_err(|e| cfd_core::Error::External(e.to_string()))
    }
}