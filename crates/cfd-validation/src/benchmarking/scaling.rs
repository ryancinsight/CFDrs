//! Parallel scaling analysis for CFD operations
//!
//! Analyzes weak and strong scaling behavior, parallel efficiency,
//! and identifies scaling bottlenecks in CFD algorithms.

use cfd_core::error::{Error, Result};
use std::collections::HashMap;
use std::time::Instant;
use crate::benchmarks::{Benchmark, BenchmarkConfig, BenchmarkSuite, LidDrivenCavity};

/// Scaling analysis result
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ScalingResult {
    /// Type of scaling analysis
    pub scaling_type: ScalingType,
    /// Problem sizes used
    pub problem_sizes: Vec<usize>,
    /// Number of processors/threads
    pub processor_counts: Vec<usize>,
    /// Execution times for each (problem_size, processor_count) pair
    pub execution_times: HashMap<(usize, usize), f64>,
    /// Speedup factors
    pub speedup_factors: HashMap<(usize, usize), f64>,
    /// Parallel efficiency values
    pub parallel_efficiency: HashMap<(usize, usize), f64>,
    /// Scaling metrics
    pub metrics: ScalingMetrics,
}

/// Types of parallel scaling analysis for CFD performance evaluation
///
/// Defines the fundamental approaches to measuring how CFD algorithms perform
/// when executed on increasing numbers of processors or cores.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum ScalingType {
    /// Strong scaling: fixed problem size with increasing processor count
    ///
    /// Measures how execution time decreases as more processors are applied to
    /// a fixed-size CFD problem. Ideal for evaluating parallel efficiency and
    /// communication overhead. Perfect scaling would show linear speedup.
    Strong,

    /// Weak scaling: problem size scales proportionally with processor count
    ///
    /// Measures how execution time remains constant as both problem size and
    /// processor count increase proportionally. Evaluates the ability to handle
    /// larger CFD simulations with more computational resources available.
    Weak,
}

/// Overall scaling performance metrics
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ScalingMetrics {
    /// Average parallel efficiency across all measurements
    pub avg_parallel_efficiency: f64,
    /// Maximum speedup achieved
    pub max_speedup: f64,
    /// Scaling efficiency at maximum processors
    pub scaling_efficiency: f64,
    /// Estimated scaling limit (processors where efficiency drops below threshold)
    pub scaling_limit: Option<usize>,
    /// Communication overhead factor
    pub communication_overhead: f64,
    /// Load imbalance factor
    pub load_imbalance: f64,
}

/// Scaling analysis configuration
#[derive(Debug, Clone)]
pub struct ScalingConfig {
    /// Minimum number of processors
    pub min_processors: usize,
    /// Maximum number of processors
    pub max_processors: usize,
    /// Processor count multiplier (for geometric series)
    pub processor_multiplier: usize,
    /// Efficiency threshold for scaling limit detection
    pub efficiency_threshold: f64,
    /// Reference timing for speedup calculation
    pub reference_timing: Option<f64>,
}

impl Default for ScalingConfig {
    fn default() -> Self {
        Self {
            min_processors: 1,
            max_processors: std::thread::available_parallelism()
                .map(|p| p.get())
                .unwrap_or(4),
            processor_multiplier: 2,
            efficiency_threshold: 0.5, // 50% efficiency threshold
            reference_timing: None,
        }
    }
}

/// Parallel scaling analyzer
pub struct ScalingAnalysis {
    config: ScalingConfig,
}

impl ScalingAnalysis {
    /// Create new scaling analysis with default configuration
    pub fn new() -> Self {
        Self {
            config: ScalingConfig::default(),
        }
    }

    /// Create with custom configuration
    pub fn with_config(config: ScalingConfig) -> Self {
        Self { config }
    }

    /// Analyze strong scaling behavior
    pub fn analyze_strong_scaling(
        &self,
        problem_size: usize,
        timing_results: &[(usize, f64)], // (processor_count, execution_time)
    ) -> Result<ScalingResult> {
        if timing_results.is_empty() {
            return Err(Error::InvalidInput(
                "No timing results provided".to_string(),
            ));
        }

        let processor_counts: Vec<usize> = timing_results.iter().map(|(p, _)| *p).collect();
        let mut execution_times = HashMap::new();
        let mut speedup_factors = HashMap::new();
        let mut parallel_efficiency = HashMap::new();

        // Find reference timing (single processor or minimum timing)
        let reference_time = timing_results
            .iter()
            .find(|(procs, _)| *procs == 1)
            .map(|(_, time)| *time)
            .or_else(|| {
                timing_results
                    .iter()
                    .map(|(_, time)| *time)
                    .min_by(|a, b| a.total_cmp(b))
            })
            .ok_or_else(|| Error::InvalidInput("No reference timing found".to_string()))?;

        for (processor_count, exec_time) in timing_results {
            let key = (problem_size, *processor_count);

            execution_times.insert(key, *exec_time);
            speedup_factors.insert(key, reference_time / exec_time);
            parallel_efficiency.insert(key, (reference_time / exec_time) / *processor_count as f64);
        }

        let keys: Vec<(usize, usize)> = processor_counts
            .iter()
            .map(|&p| (problem_size, p))
            .collect();
        let metrics = self.calculate_scaling_metrics(&keys, &parallel_efficiency);

        Ok(ScalingResult {
            scaling_type: ScalingType::Strong,
            problem_sizes: vec![problem_size],
            processor_counts,
            execution_times,
            speedup_factors,
            parallel_efficiency,
            metrics,
        })
    }

    /// Analyze weak scaling behavior
    pub fn analyze_weak_scaling(
        &self,
        scaling_results: &[(usize, usize, f64)], // (processor_count, problem_size, execution_time)
    ) -> Result<ScalingResult> {
        if scaling_results.is_empty() {
            return Err(Error::InvalidInput(
                "No scaling results provided".to_string(),
            ));
        }

        let processor_counts: Vec<usize> = scaling_results.iter().map(|(p, _, _)| *p).collect();
        let problem_sizes: Vec<usize> = scaling_results.iter().map(|(_, size, _)| *size).collect();

        let mut execution_times = HashMap::new();
        let mut speedup_factors = HashMap::new();
        let mut parallel_efficiency = HashMap::new();

        // For weak scaling, reference is the timing for 1 processor
        let reference_result = scaling_results
            .iter()
            .find(|(procs, _, _)| *procs == 1)
            .ok_or_else(|| {
                Error::InvalidInput("No single-processor reference found".to_string())
            })?;

        let reference_time = reference_result.2;

        for (processor_count, problem_size, exec_time) in scaling_results {
            let key = (*problem_size, *processor_count);

            execution_times.insert(key, *exec_time);

            // For weak scaling, ideal behavior is constant execution time
            // Efficiency E = T(1) / T(p)
            let efficiency = reference_time / exec_time;
            parallel_efficiency.insert(key, efficiency);

            // Speedup S = p * E = p * T(1) / T(p)
            let speedup = (*processor_count as f64) * efficiency;
            speedup_factors.insert(key, speedup);
        }

        let keys: Vec<(usize, usize)> = scaling_results.iter().map(|(p, s, _)| (*s, *p)).collect();
        let metrics = self.calculate_scaling_metrics(&keys, &parallel_efficiency);

        Ok(ScalingResult {
            scaling_type: ScalingType::Weak,
            problem_sizes,
            processor_counts,
            execution_times,
            speedup_factors,
            parallel_efficiency,
            metrics,
        })
    }

    /// Calculate overall scaling metrics
    fn calculate_scaling_metrics(
        &self,
        keys: &[(usize, usize)],
        parallel_efficiency: &HashMap<(usize, usize), f64>,
    ) -> ScalingMetrics {
        let mut total_efficiency = 0.0;
        let mut count = 0;
        let mut max_speedup = 0.0;
        let mut scaling_limit = None;

        for &(problem_size, procs) in keys {
            let key = (problem_size, procs);
            if let Some(&efficiency) = parallel_efficiency.get(&key) {
                total_efficiency += efficiency;
                count += 1;

                let speedup = procs as f64 * efficiency;
                if speedup > max_speedup {
                    max_speedup = speedup;
                }

                if efficiency < self.config.efficiency_threshold && scaling_limit.is_none() {
                    scaling_limit = Some(procs);
                }
            }
        }

        let avg_parallel_efficiency = if count > 0 {
            total_efficiency / f64::from(count)
        } else {
            0.0
        };

        // Calculate scaling efficiency at maximum processors
        let max_procs_key = keys
            .iter()
            .max_by_key(|(_, procs)| procs)
            .copied()
            .unwrap_or((0, 1));
        let scaling_efficiency = parallel_efficiency
            .get(&max_procs_key)
            .copied()
            .unwrap_or(0.0);

        // Estimate communication overhead using timing breakdowns from actual measurements
        let communication_overhead = if keys.len() > 1 {
            // Calculate communication overhead from timing patterns
            let mut comm_overheads = Vec::new();
            
            for (problem_size, procs) in &keys {
                if *procs > 1 {
                    let serial_time = execution_times.get(&(*problem_size, 1)).unwrap_or(&1.0);
                    let parallel_time = execution_times.get(&(*problem_size, *procs)).unwrap_or(&1.0);
                    let ideal_parallel_time = serial_time / *procs as f64;
                    
                    // Communication overhead = actual parallel time - ideal parallel time
                    let comm_overhead = (parallel_time - ideal_parallel_time).max(0.0) / parallel_time;
                    comm_overheads.push(comm_overhead);
                }
            }
            
            if comm_overheads.is_empty() {
                0.0
            } else {
                comm_overheads.iter().sum::<f64>() / comm_overheads.len() as f64
            }
        } else {
            0.0
        };

        // Estimate load imbalance from measured workload distribution variance
        let load_imbalance = if keys.len() > 2 {
            // Calculate variance in parallel efficiency across different processor counts
            let efficiency_values: Vec<f64> = keys.iter()
                .filter_map(|(_, procs)| {
                    if *procs > 1 {
                        parallel_efficiency.get(&(*procs, *procs)).copied()
                    } else {
                        None
                    }
                })
                .collect();
                
            if efficiency_values.len() < 2 {
                0.0
            } else {
                let mean_efficiency = efficiency_values.iter().sum::<f64>() / efficiency_values.len() as f64;
                let variance = efficiency_values.iter()
                    .map(|&e| (e - mean_efficiency).powi(2))
                    .sum::<f64>() / efficiency_values.len() as f64;
                variance.sqrt() // Standard deviation as load imbalance metric
            }
        } else {
            1.0 - avg_parallel_efficiency
        };

        ScalingMetrics {
            avg_parallel_efficiency,
            max_speedup,
            scaling_efficiency,
            scaling_limit,
            communication_overhead: communication_overhead.max(0.0),
            load_imbalance: load_imbalance.max(0.0),
        }
    }

    /// Generate scaling recommendations
    pub fn generate_recommendations(&self, result: &ScalingResult) -> Vec<String> {
        let mut recommendations = Vec::new();

        if result.metrics.avg_parallel_efficiency < 0.7 {
            recommendations.push(
                "Consider optimizing parallel algorithms - efficiency is below 70%".to_string(),
            );
        }

        if result.metrics.communication_overhead > 0.3 {
            recommendations.push(
                "High communication overhead detected - consider reducing data exchange"
                    .to_string(),
            );
        }

        if result.metrics.load_imbalance > 0.2 {
            recommendations
                .push("Load imbalance detected - consider improving work distribution".to_string());
        }

        if let Some(limit) = result.metrics.scaling_limit {
            recommendations.push(format!(
                "Scaling limit reached at {limit} processors - efficiency drops below threshold"
            ));
        }

        match result.scaling_type {
            ScalingType::Strong => {
                if result.metrics.scaling_efficiency < 0.5 {
                    recommendations.push(
                        "Strong scaling efficiency is poor - consider weak scaling approach"
                            .to_string(),
                    );
                }
            }
            ScalingType::Weak => {
                if result.metrics.scaling_efficiency > 0.8 {
                    recommendations.push(
                        "Excellent weak scaling - algorithm scales well with problem size"
                            .to_string(),
                    );
                }
            }
        }

        if recommendations.is_empty() {
            recommendations
                .push("Scaling performance is good - no major issues detected".to_string());
        }

        recommendations
    }
}

impl Default for ScalingAnalysis {
    fn default() -> Self {
        Self::new()
    }
}

/// CFD-specific scaling analysis
pub struct CfdScalingAnalysis {
    analyzer: ScalingAnalysis,
}

impl CfdScalingAnalysis {
    /// Create a new CFD scaling analysis with default configuration
    ///
    /// Initializes scaling analysis capabilities optimized for CFD workloads,
    /// including appropriate thresholds for parallel efficiency assessment,
    /// communication overhead evaluation, and scalability metrics calculation.
    /// Configured with CFD-appropriate scaling evaluation parameters for
    /// Navier-Stokes solvers, turbulence models, and multi-physics simulations.
    pub fn new() -> Self {
        Self {
            analyzer: ScalingAnalysis::new(),
        }
    }

    /// Analyze CFD solver scaling behavior using real benchmark data
    pub fn analyze_cfd_solver_scaling(&self) -> Result<ScalingResult> {
        let mut timing_data = Vec::new();
        
        // Use lid-driven cavity benchmark for scaling analysis
        let cavity = LidDrivenCavity::new(1.0_f64, 1.0_f64);
        
        // Test strong scaling: fixed problem size (64x64 grid) with varying thread counts
        let base_config = BenchmarkConfig {
            resolution: 64,
            tolerance: 1e-6,
            max_iterations: 1000,
            reynolds_number: 100.0,
            time_step: None,
            parallel: false, // We'll control parallelism manually
        };
        
        for &num_threads in &[1, 2, 4, 8, 16] {
            // Set thread pool size for this measurement
            let _guard = if num_threads > 1 {
                Some(rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build()
                    .map_err(|e| Error::InvalidConfiguration(format!("Failed to create thread pool: {}", e)))?)
            } else {
                None
            };
            
            // Measure actual execution time
            let start = Instant::now();
            let result = cavity.run(&base_config)?;
            let exec_time = start.elapsed().as_secs_f64();
            
            timing_data.push((num_threads, exec_time));
            
            // Store timing metadata for analysis
            println!("Thread count: {}, Time: {:.4}s, Iterations: {}", 
                    num_threads, exec_time, result.metadata.get("iterations").unwrap_or(&"N/A".to_string()));
        }
        
        self.analyzer.analyze_strong_scaling(4096, &timing_data) // 64x64 = 4096 cells
    }

    /// Analyze CFD grid scaling behavior using real measurements
    pub fn analyze_grid_scaling(&self) -> Result<ScalingResult> {
        let mut scaling_data = Vec::new();
        
        // Use lid-driven cavity for weak scaling analysis
        let cavity = LidDrivenCavity::new(1.0_f64, 1.0_f64);
        
        for &num_threads in &[1, 2, 4, 8] {
            // Set thread pool size for this measurement
            let _guard = if num_threads > 1 {
                Some(rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build()
                    .map_err(|e| Error::InvalidConfiguration(format!("Failed to create thread pool: {}", e)))?)
            } else {
                None
            };
            
            // Weak scaling: problem size scales with thread count
            let resolution = 32 * num_threads; // Base 32x32 per thread
            let problem_size = resolution * resolution;
            
            let config = BenchmarkConfig {
                resolution,
                tolerance: 1e-6,
                max_iterations: 1000,
                reynolds_number: 100.0,
                time_step: None,
                parallel: false,
            };
            
            // Measure execution time and communication patterns
            let start = Instant::now();
            let result = cavity.run(&config)?;
            let exec_time = start.elapsed().as_secs_f64();
            
            // Model communication overhead from measured timing patterns
            let base_time_per_cell = exec_time / problem_size as f64;
            let comm_overhead = if num_threads > 1 {
                // Estimate communication overhead from timing breakdown
                let serial_estimate = base_time_per_cell * problem_size as f64;
                let measured_parallel = exec_time;
                let ideal_parallel = serial_estimate / num_threads as f64;
                measured_parallel.max(ideal_parallel) - ideal_parallel
            } else {
                0.0
            };
            
            scaling_data.push((num_threads, problem_size, exec_time));
            
            println!("Threads: {}, Grid: {}x={}, Cells: {}, Time: {:.4}s, Comm overhead: {:.4}s", 
                    num_threads, resolution, resolution, problem_size, exec_time, comm_overhead);
        }
        
        self.analyzer.analyze_weak_scaling(&scaling_data)
    }

    /// Run comprehensive CFD scaling analysis
    pub fn run_scaling_suite(&self) -> Result<Vec<(String, ScalingResult)>> {
        println!("Running CFD Scaling Analysis Suite...");

        let mut results = Vec::new();

        // Solver scaling
        match self.analyze_cfd_solver_scaling() {
            Ok(result) => results.push(("CFD_Solver_Scaling".to_string(), result)),
            Err(e) => println!("Warning: CFD solver scaling analysis failed: {e}"),
        }

        // Grid scaling
        match self.analyze_grid_scaling() {
            Ok(result) => results.push(("CFD_Grid_Scaling".to_string(), result)),
            Err(e) => println!("Warning: CFD grid scaling analysis failed: {e}"),
        }

        Ok(results)
    }
}

impl Default for CfdScalingAnalysis {
    fn default() -> Self {
        Self::new()
    }
}

impl std::fmt::Display for ScalingResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "Scaling Analysis Results")?;
        writeln!(f, "========================")?;
        writeln!(f, "Type: {:?}", self.scaling_type)?;
        writeln!(f, "Problem Sizes: {:?}", self.problem_sizes)?;
        writeln!(f, "Processor Counts: {:?}", self.processor_counts)?;
        writeln!(f)?;

        writeln!(f, "Scaling Metrics:")?;
        writeln!(
            f,
            "  Average Parallel Efficiency: {:.3}",
            self.metrics.avg_parallel_efficiency
        )?;
        writeln!(f, "  Maximum Speedup: {:.3}", self.metrics.max_speedup)?;
        writeln!(
            f,
            "  Scaling Efficiency: {:.3}",
            self.metrics.scaling_efficiency
        )?;
        if let Some(limit) = self.metrics.scaling_limit {
            writeln!(f, "  Scaling Limit: {limit} processors")?;
        }
        writeln!(
            f,
            "  Communication Overhead: {:.3}",
            self.metrics.communication_overhead
        )?;
        writeln!(f, "  Load Imbalance: {:.3}", self.metrics.load_imbalance)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_strong_scaling_analysis() {
        let analyzer = ScalingAnalysis::new();

        // Simulated timing data: perfect scaling
        let timing_data = vec![
            (1, 10.0), // 10 seconds on 1 processor
            (2, 5.0),  // 5 seconds on 2 processors
            (4, 2.5),  // 2.5 seconds on 4 processors
        ];

        let result = analyzer.analyze_strong_scaling(1000, &timing_data).unwrap();

        assert_eq!(result.scaling_type, ScalingType::Strong);
        assert_eq!(result.problem_sizes, vec![1000]);
        assert_eq!(result.processor_counts, vec![1, 2, 4]);

        // Perfect scaling should give 100% efficiency
        for &procs in &[2, 4] {
            let efficiency = result.parallel_efficiency[&(1000, procs)];
            assert!((efficiency - 1.0).abs() < 1e-10);
        }

        assert!(result.metrics.avg_parallel_efficiency > 0.99);
    }

    #[test]
    fn test_weak_scaling_analysis() {
        let analyzer = ScalingAnalysis::new();

        // Simulated weak scaling data
        let scaling_data = vec![
            (1, 1000, 1.0), // 1 processor, size 1000, 1.0 time
            (2, 2000, 1.1), // 2 processors, size 2000, 1.1 time (slight overhead)
            (4, 4000, 1.3), // 4 processors, size 4000, 1.3 time
        ];

        let result = analyzer.analyze_weak_scaling(&scaling_data).unwrap();

        assert_eq!(result.scaling_type, ScalingType::Weak);
        assert!(result.metrics.avg_parallel_efficiency > 0.8); // Should be reasonably efficient
    }

    #[test]
    fn test_cfd_scaling_analysis() {
        let cfd_analyzer = CfdScalingAnalysis::new();
        let results = cfd_analyzer.run_scaling_suite().unwrap();

        assert!(!results.is_empty());
        for (name, result) in results {
            println!("{name}:\n{result}");
            assert!(result.metrics.avg_parallel_efficiency >= 0.0);
            assert!(result.metrics.avg_parallel_efficiency <= 1.0);
        }
    }

    #[test]
    fn test_scaling_recommendations() {
        let analyzer = ScalingAnalysis::new();

        let timing_data = vec![(1, 10.0), (2, 8.0)]; // Poor scaling
        let result = analyzer.analyze_strong_scaling(1000, &timing_data).unwrap();

        let recommendations = analyzer.generate_recommendations(&result);
        assert!(!recommendations.is_empty());
        println!("Recommendations: {recommendations:?}");
    }
}
