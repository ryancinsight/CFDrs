//! Supporting data structures for performance validation.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;
use std::time::Duration;

/// Communication overhead analysis results
#[derive(Debug, Clone)]
pub struct CommunicationAnalysis {
    /// Total communication time
    pub total_comm_time: Duration,
    /// Communication time per process
    pub comm_time_per_process: Vec<Duration>,
    /// Message sizes exchanged
    pub message_sizes: Vec<usize>,
    /// Network bandwidth utilization
    pub bandwidth_utilization: f64,
    /// Communication pattern efficiency
    pub pattern_efficiency: f64,
    /// Synchronization overhead
    pub sync_overhead: Duration,
    /// Load imbalance impact on communication
    pub load_imbalance_impact: f64,
}

/// Production readiness assessment
#[derive(Debug, Clone)]
pub struct ProductionReadinessReport {
    /// Overall readiness score (0-100)
    pub overall_score: u8,
    /// Component readiness scores
    pub component_scores: HashMap<String, u8>,
    /// Critical issues preventing production deployment
    pub critical_issues: Vec<String>,
    /// Performance bottlenecks
    pub bottlenecks: Vec<String>,
    /// Optimization opportunities
    pub optimization_opportunities: Vec<String>,
    /// Recommended deployment configuration
    pub recommended_config: DeploymentConfig,
    /// Scaling limits and recommendations
    pub scaling_limits: ScalingLimits,
}

/// Recommended deployment configuration
#[derive(Debug, Clone)]
pub struct DeploymentConfig {
    /// Optimal number of cores per node
    pub cores_per_node: usize,
    /// Recommended MPI implementation
    pub mpi_implementation: String,
    /// Memory requirements per core
    pub memory_per_core_mb: usize,
    /// Network requirements
    pub network_requirements: String,
    /// Recommended problem sizes
    pub recommended_problem_sizes: Vec<usize>,
}

/// Scaling limits and recommendations
#[derive(Debug, Clone)]
pub struct ScalingLimits {
    /// Maximum recommended cores
    pub max_cores: usize,
    /// Efficiency degradation threshold
    pub efficiency_threshold: f64,
    /// Communication overhead limit
    pub comm_overhead_limit: f64,
    /// Memory scaling limits
    pub memory_limits: String,
}

/// Load balancing validation results
#[derive(Debug, Clone)]
pub struct LoadBalancingValidation {
    /// Initial load imbalance ratio
    pub initial_imbalance: f64,
    /// Final load imbalance ratio after balancing
    pub final_imbalance: f64,
    /// Reduction in imbalance ratio
    pub imbalance_reduction: f64,
    /// Performance improvement ratio (balanced_time / initial_time)
    pub performance_improvement: f64,
    /// Overall effectiveness score (0-1)
    pub effectiveness_score: f64,
}

/// Simulation data for performance analysis
#[derive(Debug, Clone)]
pub struct SimulationData<T: RealField + Copy + FromPrimitive + std::fmt::LowerExp> {
    /// Total computation time
    pub comp_time: Duration,
    /// Total communication time
    pub comm_time: Duration,
    /// Total I/O time
    pub io_time: Duration,
    /// Total data transferred in bytes
    pub total_data_transferred: usize,
    /// Problem size
    pub problem_size: usize,
    /// Number of iterations
    pub iterations: usize,
    /// Final residual
    pub final_residual: T,
    /// Convergence achieved
    pub converged: bool,
}
