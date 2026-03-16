//! Performance benchmarking framework for CFD operations
//!
//! Provides comprehensive benchmarking capabilities including:
//! - Timing measurements for CFD operations
//! - Memory usage profiling
//! - Parallel scaling analysis
//! - Performance regression detection
//! - Automated benchmarking suites

pub mod analysis;
pub mod config;
pub mod memory;
pub mod performance;
pub mod production;
pub mod scaling;
pub mod suite;
pub mod timing;
pub mod visualization;

pub use analysis::{
    AlertSeverity, PerformanceAnalyzer, PerformanceReport, PerformanceTrend, RegressionAlert,
    RegressionConfig, TrendType,
};
pub use config::BenchmarkConfig;
pub use memory::{MemoryProfiler, MemoryStats};
pub use performance::{
    AlgorithmComplexity, CfdPerformanceBenchmarks, PerformanceBenchmark, PerformanceProfile,
    TimingResult,
};
pub use scaling::{ScalingAnalysis, ScalingResult};
pub use suite::{BenchmarkResult, BenchmarkStatus, BenchmarkSuite};
pub use timing::{BenchmarkStats, BenchmarkTimer};


