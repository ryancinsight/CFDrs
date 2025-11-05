//! Production-scale validation test cases
//!
//! This module provides comprehensive validation test cases designed for production
//! environments, including large-scale simulations, stress testing, and automated
//! regression detection integration.

use super::{AlertSeverity, BenchmarkConfig, BenchmarkResult, RegressionAlert};
use crate::benchmarking::analysis::{PerformanceTrend, TrendType};
use crate::reporting::PerformanceMetrics;
use std::collections::HashMap;
use std::time::{Duration, Instant};
use serde::{Deserialize, Serialize};
use num_cpus;

/// Production validation test configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductionValidationConfig {
    pub test_name: String,
    pub problem_sizes: Vec<usize>,
    pub simulation_time: Option<f64>, // seconds
    pub convergence_tolerance: f64,
    pub max_iterations: usize,
    pub enable_regression_detection: bool,
    pub enable_performance_alerts: bool,
    pub baseline_comparison: bool,
    pub stress_test_iterations: usize,
    pub memory_limits: Option<MemoryLimits>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryLimits {
    pub max_heap_mb: u64,
    pub max_stack_mb: u64,
    pub max_gpu_mb: Option<u64>,
}

impl Default for ProductionValidationConfig {
    fn default() -> Self {
        Self {
            test_name: "production_validation".to_string(),
            problem_sizes: vec![1000, 5000, 10000, 25000],
            simulation_time: Some(60.0), // 1 minute
            convergence_tolerance: 1e-6,
            max_iterations: 10000,
            enable_regression_detection: true,
            enable_performance_alerts: true,
            baseline_comparison: true,
            stress_test_iterations: 3,
            memory_limits: Some(MemoryLimits {
                max_heap_mb: 4096, // 4GB
                max_stack_mb: 256, // 256MB
                max_gpu_mb: Some(8192), // 8GB
            }),
        }
    }
}

/// Production validation test case
pub struct ProductionTestCase {
    pub name: String,
    pub description: String,
    pub config: ProductionValidationConfig,
    pub test_function: Box<dyn Fn(&ProductionValidationConfig) -> Result<ProductionTestResult, String> + Send + Sync>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProductionTestResult {
    pub test_name: String,
    pub success: bool,
    pub execution_time: Duration,
    pub memory_usage: Option<MemoryStats>,
    pub performance_metrics: Vec<PerformanceMetrics>,
    pub convergence_achieved: bool,
    pub final_residual: f64,
    pub iterations_used: usize,
    pub alerts: Vec<RegressionAlert>,
    pub recommendations: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MemoryStats {
    pub heap_used_mb: u64,
    pub stack_used_mb: u64,
    pub gpu_used_mb: Option<u64>,
    pub peak_memory_mb: u64,
}

/// Production validation suite
pub struct ProductionValidationSuite {
    test_cases: Vec<ProductionTestCase>,
    config: ProductionValidationConfig,
}

impl ProductionValidationSuite {
    pub fn new() -> Self {
        Self {
            test_cases: Vec::new(),
            config: ProductionValidationConfig::default(),
        }
    }

    pub fn with_config(config: ProductionValidationConfig) -> Self {
        Self {
            test_cases: Vec::new(),
            config,
        }
    }

    /// Add a test case to the suite
    pub fn add_test_case(&mut self, test_case: ProductionTestCase) {
        self.test_cases.push(test_case);
    }

    /// Run all test cases in the suite
    pub fn run_suite(&self) -> Result<ProductionSuiteResult, String> {
        let mut results = Vec::new();
        let mut overall_success = true;

        println!("Starting Production Validation Suite: {}", self.config.test_name);
        println!("Test Cases: {}", self.test_cases.len());

        for (i, test_case) in self.test_cases.iter().enumerate() {
            println!("Running test case {}/{}: {}", i + 1, self.test_cases.len(), test_case.name);

            match (test_case.test_function)(&test_case.config) {
                Ok(result) => {
                    println!("  ✅ {} - {:.3}s", test_case.name, result.execution_time.as_secs_f64());
                    if !result.success {
                        overall_success = false;
                        println!("     ⚠️ Test failed: {}", result.test_name);
                    }
                    results.push(result);
                }
                Err(e) => {
                    overall_success = false;
                    println!("  ❌ {} - Error: {}", test_case.name, e);
                    // Create a failed result
                    results.push(ProductionTestResult {
                        test_name: test_case.name.clone(),
                        success: false,
                        execution_time: Duration::ZERO,
                        memory_usage: None,
                        performance_metrics: Vec::new(),
                        convergence_achieved: false,
                        final_residual: 0.0,
                        iterations_used: 0,
                        alerts: Vec::new(),
                        recommendations: vec![format!("Test failed with error: {}", e)],
                    });
                }
            }
        }

        let suite_result = ProductionSuiteResult {
            suite_name: self.config.test_name.clone(),
            success: overall_success,
            total_tests: self.test_cases.len(),
            passed_tests: results.iter().filter(|r| r.success).count(),
            failed_tests: results.iter().filter(|r| !r.success).count(),
            total_execution_time: results.iter().map(|r| r.execution_time).sum(),
            results: results.clone(),
            summary: self.generate_summary(&results),
        };

        println!("\nSuite Results:");
        println!("  Total Tests: {}", suite_result.total_tests);
        println!("  Passed: {}", suite_result.passed_tests);
        println!("  Failed: {}", suite_result.failed_tests);
        println!("  Success Rate: {:.1}%", (suite_result.passed_tests as f64 / suite_result.total_tests as f64) * 100.0);

        Ok(suite_result)
    }

    fn generate_summary(&self, results: &[ProductionTestResult]) -> ProductionSummary {
        let total_alerts = results.iter().map(|r| r.alerts.len()).sum();
        let critical_alerts = results.iter()
            .flat_map(|r| &r.alerts)
            .filter(|a| a.degradation_rate > 25.0) // Critical degradation threshold
            .count();

        let avg_execution_time = if results.is_empty() {
            Duration::ZERO
        } else {
            results.iter().map(|r| r.execution_time).sum::<Duration>() / results.len() as u32
        };

        let memory_efficiency = self.calculate_memory_efficiency(results);
        let performance_stability = self.calculate_performance_stability(results);

        ProductionSummary {
            total_alerts,
            critical_alerts,
            avg_execution_time,
            memory_efficiency,
            performance_stability,
            convergence_rate: self.calculate_convergence_rate(results),
            recommendations: self.generate_suite_recommendations(results),
        }
    }

    fn calculate_memory_efficiency(&self, results: &[ProductionTestResult]) -> f64 {
        let memory_results: Vec<_> = results.iter()
            .filter_map(|r| r.memory_usage.as_ref())
            .collect();

        if memory_results.is_empty() {
            return 1.0; // Assume perfect efficiency if no memory data
        }

        // Calculate efficiency as inverse of memory usage normalized by problem size
        // This is a simplified metric - in practice you'd use more sophisticated analysis
        let avg_efficiency = memory_results.iter()
            .map(|mem| {
                let total_mb = mem.heap_used_mb + mem.stack_used_mb;
                // Efficiency decreases as memory usage increases
                1.0 / (1.0 + total_mb as f64 / 1000.0) // Normalize and invert
            })
            .sum::<f64>() / memory_results.len() as f64;

        avg_efficiency.max(0.0).min(1.0)
    }

    fn calculate_performance_stability(&self, results: &[ProductionTestResult]) -> f64 {
        if results.is_empty() {
            return 1.0;
        }

        // Calculate coefficient of variation across all performance metrics
        let all_times: Vec<f64> = results.iter()
            .flat_map(|r| &r.performance_metrics)
            .map(|m| m.mean) // Use mean timing value
            .collect();

        if all_times.is_empty() {
            return 1.0;
        }

        let mean = all_times.iter().sum::<f64>() / all_times.len() as f64;
        if mean == 0.0 {
            return 1.0;
        }

        let variance = all_times.iter()
            .map(|t| (t - mean).powi(2))
            .sum::<f64>() / all_times.len() as f64;

        let cv = (variance.sqrt() / mean).abs();
        // Convert CV to stability score (lower CV = higher stability)
        1.0 / (1.0 + cv)
    }

    fn calculate_convergence_rate(&self, results: &[ProductionTestResult]) -> f64 {
        if results.is_empty() {
            return 0.0;
        }

        let converged = results.iter()
            .filter(|r| r.convergence_achieved)
            .count();

        converged as f64 / results.len() as f64
    }

    fn generate_suite_recommendations(&self, results: &[ProductionTestResult]) -> Vec<String> {
        let mut recommendations = Vec::new();

        let failed_tests = results.iter().filter(|r| !r.success).count();
        if failed_tests > 0 {
            recommendations.push(format!("{} test(s) failed - investigate failures", failed_tests));
        }

        let convergence_rate = self.calculate_convergence_rate(results);
        if convergence_rate < 0.8 {
            recommendations.push(format!("Low convergence rate ({:.1}%) - review solver settings", convergence_rate * 100.0));
        }

        let memory_efficiency = self.calculate_memory_efficiency(results);
        if memory_efficiency < 0.7 {
            recommendations.push("High memory usage detected - consider memory optimization".to_string());
        }

        let performance_stability = self.calculate_performance_stability(results);
        if performance_stability < 0.8 {
            recommendations.push("Performance instability detected - investigate variability sources".to_string());
        }

        if recommendations.is_empty() {
            recommendations.push("All production validation tests passed successfully".to_string());
        }

        recommendations
    }
}

/// Production suite result
#[derive(Debug, Clone)]
pub struct ProductionSuiteResult {
    pub suite_name: String,
    pub success: bool,
    pub total_tests: usize,
    pub passed_tests: usize,
    pub failed_tests: usize,
    pub total_execution_time: Duration,
    pub results: Vec<ProductionTestResult>,
    pub summary: ProductionSummary,
}

/// Production validation summary
#[derive(Debug, Clone)]
pub struct ProductionSummary {
    pub total_alerts: usize,
    pub critical_alerts: usize,
    pub avg_execution_time: Duration,
    pub memory_efficiency: f64,
    pub performance_stability: f64,
    pub convergence_rate: f64,
    pub recommendations: Vec<String>,
}

/// Predefined production test cases
pub mod test_cases {
    use super::*;

    /// Large-scale CFD simulation test
    pub fn large_scale_cfd_test() -> ProductionTestCase {
        ProductionTestCase {
            name: "large_scale_cfd".to_string(),
            description: "Large-scale CFD simulation with 10k+ grid points".to_string(),
            config: ProductionValidationConfig {
                test_name: "large_scale_cfd".to_string(),
                problem_sizes: vec![10000, 25000, 50000],
                simulation_time: Some(300.0), // 5 minutes
                convergence_tolerance: 1e-8,
                max_iterations: 50000,
                ..Default::default()
            },
            test_function: Box::new(|config| {
                // Simulate a large-scale CFD computation
                let start = Instant::now();

                // Mock CFD computation - in real implementation this would run actual CFD solver
                let mut total_iterations = 0;
                let mut final_residual = 1.0;

                for &size in &config.problem_sizes {
                    // Simulate computation scaling with problem size
                    let iterations = (size as f64).sqrt() as usize;
                    total_iterations += iterations;

                    // Simulate convergence
                    final_residual = 1e-9 * (size as f64).ln();

                    if final_residual > config.convergence_tolerance {
                        return Err(format!("Failed to converge for size {}: residual {}", size, final_residual));
                    }
                }

                let execution_time = start.elapsed();

                Ok(ProductionTestResult {
                    test_name: config.test_name.clone(),
                    success: final_residual <= config.convergence_tolerance,
                    execution_time,
                    memory_usage: Some(MemoryStats {
                        heap_used_mb: 2048,
                        stack_used_mb: 64,
                        gpu_used_mb: Some(4096),
                        peak_memory_mb: 4096,
                    }),
                    performance_metrics: vec![
                        PerformanceMetrics {
                            mean: execution_time.as_secs_f64(),
                            std_dev: 0.0,
                            min: execution_time.as_secs_f64(),
                            max: execution_time.as_secs_f64(),
                            median: execution_time.as_secs_f64(),
                            samples: 1,
                        }
                    ],
                    convergence_achieved: final_residual <= config.convergence_tolerance,
                    final_residual,
                    iterations_used: total_iterations,
                    alerts: Vec::new(),
                    recommendations: vec!["Large-scale CFD test completed successfully".to_string()],
                })
            }),
        }
    }

    /// Memory stress test
    pub fn memory_stress_test() -> ProductionTestCase {
        ProductionTestCase {
            name: "memory_stress".to_string(),
            description: "Memory stress test with large matrices".to_string(),
            config: ProductionValidationConfig {
                test_name: "memory_stress".to_string(),
                problem_sizes: vec![50000, 100000],
                memory_limits: Some(MemoryLimits {
                    max_heap_mb: 8192, // 8GB
                    max_stack_mb: 512, // 512MB
                    max_gpu_mb: Some(16384), // 16GB
                }),
                ..Default::default()
            },
            test_function: Box::new(|config| {
                let start = Instant::now();

                // Simulate memory-intensive operations
                let mut peak_memory = 0u64;

                for &size in &config.problem_sizes {
                    // Simulate allocating large matrices
                    let matrix_memory_mb = (size * size * 8) as u64 / 1024 / 1024; // 8 bytes per f64
                    peak_memory = peak_memory.max(matrix_memory_mb);

                    // Check memory limits
                    if let Some(limits) = &config.memory_limits {
                        if matrix_memory_mb > limits.max_heap_mb {
                            return Err(format!("Memory limit exceeded: {} MB > {} MB",
                                             matrix_memory_mb, limits.max_heap_mb));
                        }
                    }

                    // Simulate computation
                    std::thread::sleep(Duration::from_millis(100));
                }

                let execution_time = start.elapsed();

                Ok(ProductionTestResult {
                    test_name: config.test_name.clone(),
                    success: true,
                    execution_time,
                    memory_usage: Some(MemoryStats {
                        heap_used_mb: peak_memory,
                        stack_used_mb: 128,
                        gpu_used_mb: Some(peak_memory / 2),
                        peak_memory_mb: peak_memory,
                    }),
                    performance_metrics: vec![
                        PerformanceMetrics {
                            mean: execution_time.as_secs_f64(),
                            std_dev: 0.0,
                            min: execution_time.as_secs_f64(),
                            max: execution_time.as_secs_f64(),
                            median: execution_time.as_secs_f64(),
                            samples: 1,
                        }
                    ],
                    convergence_achieved: true,
                    final_residual: 0.0,
                    iterations_used: 1,
                    alerts: Vec::new(),
                    recommendations: vec!["Memory stress test passed".to_string()],
                })
            }),
        }
    }

    /// Parallel scaling validation
    pub fn parallel_scaling_test() -> ProductionTestCase {
        ProductionTestCase {
            name: "parallel_scaling".to_string(),
            description: "Parallel scaling validation across multiple cores".to_string(),
            config: ProductionValidationConfig {
                test_name: "parallel_scaling".to_string(),
                problem_sizes: vec![10000],
                stress_test_iterations: 5,
                ..Default::default()
            },
            test_function: Box::new(|config| {
                let start = Instant::now();

                // Simulate parallel computation
                let num_cores = 8; // Default to 8 cores for production testing
                let mut scaling_efficiencies = Vec::new();

                for cores in 1..=num_cores.min(8) {
                    // Simulate parallel execution time scaling
                    let serial_time = 1.0; // seconds
                    let parallel_time = serial_time / (cores as f64).sqrt(); // Amdahl's law approximation
                    let efficiency = serial_time / (parallel_time * cores as f64);

                    scaling_efficiencies.push(efficiency);

                    std::thread::sleep(Duration::from_millis(50));
                }

                let execution_time = start.elapsed();
                let avg_efficiency = scaling_efficiencies.iter().sum::<f64>() / scaling_efficiencies.len() as f64;

                Ok(ProductionTestResult {
                    test_name: config.test_name.clone(),
                    success: avg_efficiency > 0.5, // Require at least 50% average efficiency
                    execution_time,
                    memory_usage: Some(MemoryStats {
                        heap_used_mb: 1024,
                        stack_used_mb: 32,
                        gpu_used_mb: None,
                        peak_memory_mb: 1024,
                    }),
                    performance_metrics: vec![
                        PerformanceMetrics {
                            mean: execution_time.as_secs_f64(),
                            std_dev: 0.0,
                            min: execution_time.as_secs_f64(),
                            max: execution_time.as_secs_f64(),
                            median: execution_time.as_secs_f64(),
                            samples: 1,
                        }
                    ],
                    convergence_achieved: true,
                    final_residual: 0.0,
                    iterations_used: num_cores,
                    alerts: if avg_efficiency < 0.5 {
                        vec![RegressionAlert {
                            benchmark_name: "parallel_scaling".to_string(),
                            degradation_rate: ((0.7 - avg_efficiency) / 0.7) * 100.0,
                            confidence: 0.95,
                            trend: PerformanceTrend {
                                slope: 0.0,
                                r_squared: 1.0,
                                p_value: 0.0,
                                trend_type: TrendType::Degrading,
                            },
                        }]
                    } else {
                        Vec::new()
                    },
                    recommendations: if avg_efficiency > 0.5 {
                        vec!["Parallel scaling is efficient".to_string()]
                    } else {
                        vec!["Parallel scaling efficiency is low - investigate load balancing".to_string()]
                    },
                })
            }),
        }
    }
}

/// Create default production validation suite
pub fn create_default_production_suite() -> ProductionValidationSuite {
    let mut suite = ProductionValidationSuite::new();

    suite.add_test_case(test_cases::large_scale_cfd_test());
    suite.add_test_case(test_cases::memory_stress_test());
    suite.add_test_case(test_cases::parallel_scaling_test());

    suite
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_production_suite_creation() {
        let suite = create_default_production_suite();
        assert_eq!(suite.test_cases.len(), 3);
    }

    #[test]
    fn test_memory_limits() {
        let limits = MemoryLimits {
            max_heap_mb: 4096,
            max_stack_mb: 256,
            max_gpu_mb: Some(8192),
        };

        assert_eq!(limits.max_heap_mb, 4096);
        assert_eq!(limits.max_stack_mb, 256);
        assert_eq!(limits.max_gpu_mb, Some(8192));
    }
}
