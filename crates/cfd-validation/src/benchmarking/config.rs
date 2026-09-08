//! Configuration for CFD performance benchmarks
//!
//! Domain: CFD Performance Benchmarking
//! Bounded Context: Configuration
//! Architectural Pattern: Value Object

use std::collections::HashMap;

/// Benchmark configuration with architectural purity
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BenchmarkConfig {
    /// Number of iterations for each benchmark
    pub iterations: usize,
    /// Number of warmup iterations
    pub warmup_iterations: usize,
    /// Enable memory profiling
    pub enable_memory: bool,
    /// Enable scaling analysis
    pub enable_scaling: bool,
    /// Enable detailed reporting
    pub detailed_reporting: bool,
    /// Performance regression threshold (percentage)
    pub regression_threshold: f64,
    /// Baseline performance data (for regression detection)
    pub baseline_data: Option<HashMap<String, f64>>,
    /// Problem sizes to benchmark
    pub problem_sizes: Vec<usize>,
    /// Output directory for reports
    pub output_dir: Option<String>,
}

impl Default for BenchmarkConfig {
    fn default() -> Self {
        Self {
            iterations: 5,
            warmup_iterations: 1,
            enable_memory: true,
            enable_scaling: true,
            detailed_reporting: true,
            regression_threshold: 5.0, // 5% regression threshold
            baseline_data: None,
            problem_sizes: vec![1000, 10000, 100000],
            output_dir: Some("performance_results".to_string()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_default() {
        let config = BenchmarkConfig::default();
        assert_eq!(config.iterations, 5);
        assert!(config.enable_memory);
        assert!(config.enable_scaling);
        assert_eq!(config.problem_sizes, vec![1000, 10000, 100000]);
    }
}
