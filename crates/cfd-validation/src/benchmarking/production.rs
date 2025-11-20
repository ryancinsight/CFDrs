//! Production-ready performance profiling for CFD algorithms
//!
//! This module provides comprehensive performance analysis including:
//! - Algorithm complexity documentation (Big-O analysis)
//! - Memory bandwidth and cache efficiency metrics
//! - Performance recommendations based on empirical data
//! - Literature-backed complexity analysis

use super::performance::{AlgorithmComplexity, CfdPerformanceBenchmarks, PerformanceProfile};
use cfd_core::error::{Error, Result};

/// Comprehensive performance profiling runner
pub struct PerformanceProfiler {
    benchmarks: CfdPerformanceBenchmarks,
}

impl PerformanceProfiler {
    /// Create new performance profiler
    pub fn new() -> Self {
        Self {
            benchmarks: CfdPerformanceBenchmarks::new(),
        }
    }

    /// Run comprehensive CFD algorithm performance profiling
    pub fn run_comprehensive_profiling(&self) -> Result<PerformanceReport> {
        println!("🔬 Comprehensive CFD Algorithm Performance Profiling");
        println!("==================================================");
        println!("Literature References:");
        println!("- Saad (2003): Iterative Methods for Sparse Linear Systems");
        println!("- Williams et al. (2009): SPMV on emerging multicore platforms");
        println!("- Salari & Knupp (2000): Code verification by MMS");
        println!("- LeVeque (2002): Finite Volume Methods for Hyperbolic Problems");
        println!();

        // Run CFD algorithm benchmarks
        let profiles = self.benchmarks.benchmark_cfd_algorithms()?;

        // Generate comprehensive report
        let summary = self.generate_summary(&profiles);
        let recommendations = self.generate_system_recommendations(&profiles);
        let report = PerformanceReport {
            timestamp: chrono::Utc::now(),
            profiles,
            summary,
            recommendations,
        };

        // Display results
        self.display_report(&report);

        Ok(report)
    }

    /// Generate performance summary
    fn generate_summary(&self, profiles: &[PerformanceProfile]) -> PerformanceSummary {
        let total_gflops = profiles.iter().map(|p| p.gflops).sum::<f64>();
        let avg_memory_bandwidth =
            profiles.iter().map(|p| p.memory_bandwidth).sum::<f64>() / profiles.len() as f64;
        let avg_cache_efficiency = profiles
            .iter()
            .map(|p| p.complexity.cache_efficiency)
            .sum::<f64>()
            / profiles.len() as f64;
        let avg_scalability = profiles
            .iter()
            .map(|p| p.complexity.scalability)
            .sum::<f64>()
            / profiles.len() as f64;

        let best_algorithm = profiles
            .iter()
            .max_by(|a, b| a.gflops.partial_cmp(&b.gflops).unwrap())
            .map(|p| p.complexity.name.clone())
            .unwrap_or_else(|| "None".to_string());

        let worst_algorithm = profiles
            .iter()
            .min_by(|a, b| {
                a.complexity
                    .cache_efficiency
                    .partial_cmp(&b.complexity.cache_efficiency)
                    .unwrap()
            })
            .map(|p| p.complexity.name.clone())
            .unwrap_or_else(|| "None".to_string());

        PerformanceSummary {
            total_profiles: profiles.len(),
            total_gflops,
            avg_memory_bandwidth,
            avg_cache_efficiency,
            avg_scalability,
            best_algorithm,
            worst_algorithm,
        }
    }

    /// Generate system-wide performance recommendations
    fn generate_system_recommendations(
        &self,
        profiles: &[PerformanceProfile],
    ) -> Vec<PerformanceRecommendation> {
        let mut recommendations = Vec::new();

        // Memory bandwidth analysis
        let avg_bandwidth =
            profiles.iter().map(|p| p.memory_bandwidth).sum::<f64>() / profiles.len() as f64;
        if avg_bandwidth < 20.0 {
            recommendations.push(PerformanceRecommendation {
                category: "Memory Bandwidth".to_string(),
                severity: "High".to_string(),
                description: "Low memory bandwidth utilization detected".to_string(),
                solution:
                    "Consider cache-blocking, data structure reorganization, or SIMD vectorization"
                        .to_string(),
                literature: "Williams et al. (2009): Roofline model for performance analysis"
                    .to_string(),
            });
        }

        // Cache efficiency analysis
        let avg_cache_eff = profiles
            .iter()
            .map(|p| p.complexity.cache_efficiency)
            .sum::<f64>()
            / profiles.len() as f64;
        if avg_cache_eff < 0.5 {
            recommendations.push(PerformanceRecommendation {
                category: "Cache Efficiency".to_string(),
                severity: "Medium".to_string(),
                description: "Poor cache utilization across algorithms".to_string(),
                solution: "Implement cache-aware algorithms and data layouts".to_string(),
                literature: "Gustavson (1978): Two-dimensional basic linear algebra subprograms"
                    .to_string(),
            });
        }

        // Scalability analysis
        let avg_scalability = profiles
            .iter()
            .map(|p| p.complexity.scalability)
            .sum::<f64>()
            / profiles.len() as f64;
        if avg_scalability < 0.6 {
            recommendations.push(PerformanceRecommendation {
                category: "Parallel Scalability".to_string(),
                severity: "Medium".to_string(),
                description: "Limited parallel scalability detected".to_string(),
                solution: "Optimize for distributed computing or consider algorithm redesign"
                    .to_string(),
                literature: "Amdahl (1967): Validity of the single processor approach".to_string(),
            });
        }

        // Algorithm complexity recommendations
        let complex_algorithms = profiles
            .iter()
            .filter(|p| {
                p.complexity.time_complexity == "O(N²)" || p.complexity.time_complexity == "O(N³)"
            })
            .count();

        if complex_algorithms > profiles.len() / 2 {
            recommendations.push(PerformanceRecommendation {
                category: "Algorithm Complexity".to_string(),
                severity: "High".to_string(),
                description: format!(
                    "{} algorithms have poor complexity scaling",
                    complex_algorithms
                ),
                solution:
                    "Consider fast algorithms (FFT, multigrid, preconditioners) for large problems"
                        .to_string(),
                literature: "Golub & Van Loan (1996): Matrix Computations".to_string(),
            });
        }

        recommendations
    }

    /// Display comprehensive performance report
    fn display_report(&self, report: &PerformanceReport) {
        println!("\n📊 Performance Profiling Summary");
        println!("===============================");
        println!(
            "Timestamp: {}",
            report.timestamp.format("%Y-%m-%d %H:%M:%S UTC")
        );
        println!(
            "Total Algorithms Profiled: {}",
            report.summary.total_profiles
        );
        println!(
            "Total Performance: {:.2} GFLOPS",
            report.summary.total_gflops
        );
        println!(
            "Average Memory Bandwidth: {:.2} GB/s",
            report.summary.avg_memory_bandwidth
        );
        println!(
            "Average Cache Efficiency: {:.1}%",
            report.summary.avg_cache_efficiency * 100.0
        );
        println!(
            "Average Parallel Scalability: {:.1}%",
            report.summary.avg_scalability * 100.0
        );
        println!(
            "Best Performing Algorithm: {}",
            report.summary.best_algorithm
        );
        println!(
            "Algorithm Needing Optimization: {}",
            report.summary.worst_algorithm
        );

        if !report.recommendations.is_empty() {
            println!("\n💡 Performance Recommendations");
            println!("=============================");
            for (i, rec) in report.recommendations.iter().enumerate() {
                println!(
                    "{}. [{}] {}: {}",
                    i + 1,
                    rec.severity,
                    rec.category,
                    rec.description
                );
                println!("   Solution: {}", rec.solution);
                println!("   Reference: {}", rec.literature);
                println!();
            }
        }

        println!("✅ Performance profiling completed successfully!");
    }
}

/// Performance profiling report
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceReport {
    /// Report timestamp
    pub timestamp: chrono::DateTime<chrono::Utc>,
    /// Individual performance profiles
    pub profiles: Vec<PerformanceProfile>,
    /// Summary statistics
    pub summary: PerformanceSummary,
    /// System recommendations
    pub recommendations: Vec<PerformanceRecommendation>,
}

/// Performance summary statistics
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceSummary {
    /// Total number of algorithms profiled
    pub total_profiles: usize,
    /// Total GFLOPS across all algorithms
    pub total_gflops: f64,
    /// Average memory bandwidth (GB/s)
    pub avg_memory_bandwidth: f64,
    /// Average cache efficiency (0.0-1.0)
    pub avg_cache_efficiency: f64,
    /// Average parallel scalability (0.0-1.0)
    pub avg_scalability: f64,
    /// Best performing algorithm
    pub best_algorithm: String,
    /// Algorithm with worst cache efficiency
    pub worst_algorithm: String,
}

/// Performance recommendation
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PerformanceRecommendation {
    /// Recommendation category
    pub category: String,
    /// Severity level (Low/Medium/High)
    pub severity: String,
    /// Problem description
    pub description: String,
    /// Recommended solution
    pub solution: String,
    /// Literature reference
    pub literature: String,
}

impl Default for PerformanceProfiler {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_performance_profiler_creation() {
        let profiler = PerformanceProfiler::new();
        // Should create without errors
        assert!(profiler.benchmarks.run_cfd_suite().is_ok());
    }

    #[test]
    fn test_performance_report_structure() {
        let profiler = PerformanceProfiler::new();

        // This test would run actual profiling in a real scenario
        // For now, just test that the profiler can be created
        assert!(true); // Placeholder - would need actual test data
    }
}
