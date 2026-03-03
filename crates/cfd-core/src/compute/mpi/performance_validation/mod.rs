//! Performance validation and scaling benchmarks for MPI parallelization.
//!
//! Provides comprehensive performance validation infrastructure for assessing
//! the efficiency and scalability of MPI-parallel CFD simulations.

/// Scaling assessment and grading logic.
mod assessment;
/// Performance metrics, scaling test results, and grading types.
mod metrics;
/// Timing utilities for performance measurement.
mod timer;
/// Supporting data structures for performance validation.
mod types;
/// MPI performance validator for scaling benchmarks and production readiness.
mod validator;

pub use assessment::ScalingAssessment;
pub use metrics::{PerformanceMetrics, ScalingGrade, ScalingTestResult, ScalingTestType};
pub use timer::PerformanceTimer;
pub use types::{
    CommunicationAnalysis, DeploymentConfig, LoadBalancingValidation, ProductionReadinessReport,
    ScalingLimits, SimulationData,
};
pub use validator::PerformanceValidator;

#[cfg(test)]
mod tests {
    use super::*;
    use std::time::Duration;

    #[test]
    fn test_performance_metrics_creation() {
        let mut metrics = PerformanceMetrics::new();
        metrics.num_processes = 4;
        metrics.problem_size = 100000;
        metrics.total_time = Duration::from_secs(10);

        assert_eq!(metrics.num_processes, 4);
        assert_eq!(metrics.problem_size, 100000);
        assert_eq!(metrics.total_time, Duration::from_secs(10));
    }

    #[test]
    fn test_performance_metrics_derived_calculation() {
        let mut metrics = PerformanceMetrics::new();
        metrics.num_processes = 4;
        metrics.total_time = Duration::from_secs(100);
        metrics.communication_time = Duration::from_secs(10);
        metrics.computation_time = Duration::from_secs(85);
        metrics.load_balancing_time = Duration::from_secs(5);

        metrics.calculate_derived_metrics();

        assert_eq!(metrics.communication_overhead_percent, 10.0);
        assert_eq!(metrics.parallel_efficiency, 0.25);
    }

    #[test]
    fn test_scaling_assessment() {
        let mut results = ScalingTestResult {
            core_counts: vec![1, 2, 4, 8, 16, 32, 64],
            metrics: Vec::new(),
            scaling_efficiency: vec![1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.82],
            communication_trend: vec![0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0],
            load_imbalance_trend: vec![1.0, 1.01, 1.02, 1.03, 1.05, 1.08, 1.15],
            test_type: ScalingTestType::Strong,
            assessment: ScalingAssessment {
                grade: ScalingGrade::F,
                efficiency_at_64_cores: 0.0,
                comm_overhead_at_64_cores: 0.0,
                load_imbalance_at_64_cores: 0.0,
                notes: Vec::new(),
                recommendations: Vec::new(),
            },
        };

        results.assess_scaling();

        assert_eq!(results.assessment.grade, ScalingGrade::A);
        assert_eq!(results.assessment.efficiency_at_64_cores, 0.82);
        assert_eq!(results.assessment.comm_overhead_at_64_cores, 12.0);
        assert_eq!(results.assessment.load_imbalance_at_64_cores, 1.15);
    }

    #[test]
    fn test_performance_timer() {
        let mut timer = PerformanceTimer::new();

        timer.start();
        std::thread::sleep(Duration::from_millis(10));
        timer.stop();

        assert!(timer.total_time() >= Duration::from_millis(10));
    }
}
