//! Data structures for Richardson extrapolation and validation

use nalgebra::RealField;

/// Result of Richardson extrapolation error analysis
#[derive(Debug, Clone)]
pub struct RichardsonResult<T: RealField + Copy> {
    /// Extrapolated solution (estimate of exact solution)
    pub extrapolated_solution: T,
    /// Estimated convergence order
    pub estimated_order: T,
    /// Discretization errors for each grid level
    pub grid_errors: Vec<T>,
    /// Observed convergence rates between grid levels
    pub convergence_rates: Vec<T>,
    /// Grid sizes used in the analysis
    pub grid_sizes: Vec<usize>,
}

/// Comprehensive boundary condition validation result
#[derive(Debug, Clone)]
pub struct BoundaryValidationResult<T: RealField + Copy> {
    /// Maximum boundary condition error
    pub max_bc_error: T,
    /// Flux continuity errors at boundaries
    pub flux_continuity_errors: Vec<T>,
    /// Compatibility check results
    pub compatibility_passed: bool,
    /// Physical consistency check results
    pub physical_consistency_passed: bool,
    /// Boundary condition types validated
    pub validated_boundaries: Vec<String>,
}

/// Result of Richardson extrapolation applied to MMS
#[derive(Debug, Clone)]
pub struct RichardsonMmsResult<T: RealField + Copy> {
    /// Grid sizes used in the study
    pub grid_sizes: Vec<T>,
    /// L2 errors for each grid
    pub l2_errors: Vec<T>,
    /// Convergence study analysis
    pub convergence_study: crate::convergence::ConvergenceStudy<T>,
    /// Richardson extrapolation results
    pub richardson_results: Vec<(T, T)>, // (extrapolated_value, estimated_order)
    /// Grid convergence indices
    pub gci_values: Vec<T>,
    /// Asymptotic range indicators
    pub is_asymptotic: Vec<bool>,
}

impl<T: RealField + Copy> RichardsonMmsResult<T> {
    /// Check if all grids are in asymptotic range
    pub fn all_asymptotic(&self) -> bool {
        self.is_asymptotic.iter().all(|&x| x)
    }

    /// Get the final extrapolated solution (most accurate estimate)
    pub fn final_extrapolated_solution(&self) -> Option<T> {
        self.richardson_results.last().map(|(val, _)| *val)
    }

    /// Get the final estimated order of accuracy
    pub fn final_estimated_order(&self) -> Option<T> {
        self.richardson_results.last().map(|(_, order)| *order)
    }
}

/// Performance profiling results
#[derive(Debug, Clone)]
pub struct PerformanceProfile {
    /// Big-O complexity documentation for algorithms
    pub algorithm_complexity: AlgorithmComplexity,
    /// Memory bandwidth analysis and access patterns
    pub memory_bandwidth_analysis: MemoryBandwidthAnalysis,
    /// Cache efficiency metrics
    pub cache_efficiency_metrics: CacheEfficiencyMetrics,
    /// Parallel scalability analysis
    pub scalability_analysis: ScalabilityAnalysis,
}

/// Algorithm complexity documentation
#[derive(Debug, Clone)]
pub struct AlgorithmComplexity {
    /// Richardson extrapolation complexity
    pub richardson_extrapolation: String,
    /// Boundary validation complexity
    pub boundary_validation: String,
    /// Convergence analysis complexity
    pub convergence_analysis: String,
    /// Error estimation complexity
    pub error_estimation: String,
    /// Manufactured solution evaluation complexity
    pub manufactured_solution_evaluation: String,
    /// Grid convergence index complexity
    pub grid_convergence_index: String,
}

/// Memory bandwidth analysis
#[derive(Debug, Clone)]
pub struct MemoryBandwidthAnalysis {
    /// Memory access patterns
    pub memory_access_patterns: String,
    /// Cache line utilization
    pub cache_line_utilization: String,
    /// Memory bandwidth requirements
    pub memory_bandwidth_requirements: String,
}

/// Cache efficiency metrics
#[derive(Debug, Clone)]
pub struct CacheEfficiencyMetrics {
    /// Spatial locality
    pub spatial_locality: String,
    /// Temporal locality
    pub temporal_locality: String,
    /// Cache miss rate
    pub cache_miss_rate: String,
}

/// Parallel scalability analysis
#[derive(Debug, Clone)]
pub struct ScalabilityAnalysis {
    /// Parallel efficiency
    pub parallel_efficiency: String,
    /// Communication overhead
    pub communication_overhead: String,
    /// Load balancing
    pub load_balancing: String,
}

/// Stability region analysis for a numerical scheme
#[derive(Debug, Clone)]
pub struct StabilityRegion<T: RealField + Copy> {
    /// Maximum CFL number for stability
    pub cfl_max: T,
    /// CFL range start
    pub cfl_range_start: T,
    /// CFL range end
    pub cfl_range_end: T,
    /// Von Neumann stability analysis results
    pub von_neumann_analysis: Option<VonNeumannAnalysis<T>>,
}

/// Von Neumann stability analysis
#[derive(Debug, Clone)]
pub struct VonNeumannAnalysis<T: RealField + Copy> {
    /// Amplification factors for different wavenumbers
    pub amplification_factors: Vec<T>,
    /// Wavenumber range analyzed
    pub wavenumber_range: (T, T),
    /// Maximum amplification factor
    pub max_amplification: T,
}

/// Numerical stability analysis results
#[derive(Debug, Clone)]
pub struct NumericalStabilityAnalysis<T: RealField + Copy> {
    /// Stability regions for different schemes
    pub stability_regions: Vec<StabilityRegion<T>>,
    /// CFL conditions for various problems
    pub cfl_conditions: Vec<String>,
    /// Dispersion analysis results
    pub dispersion_analysis: Vec<String>,
}

/// Conservation property errors
#[derive(Debug, Clone)]
pub struct ConservationErrors<T: RealField + Copy> {
    /// Mass conservation error
    pub mass_conservation_error: T,
    /// Momentum conservation error
    pub momentum_conservation_error: T,
    /// Energy conservation error
    pub energy_conservation_error: T,
    /// Angular momentum conservation error
    pub angular_momentum_error: T,
}

/// Conservation property verification
#[derive(Debug, Clone)]
pub struct ConservationAnalysis<T: RealField + Copy> {
    /// Conservation errors measured
    pub conservation_errors: Vec<ConservationErrors<T>>,
    /// Conservation properties verified
    pub conservation_properties: Vec<String>,
}

/// Edge case validation result
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckStatus {
    /// Validation passed
    Passed,
    /// Validation failed
    Failed,
}

impl CheckStatus {
    #[allow(missing_docs)]
    pub fn is_passed(&self) -> bool {
        matches!(self, Self::Passed)
    }
}

#[allow(missing_docs)]
#[derive(Debug, Clone)]
pub struct EdgeCaseResult {
    /// Boundary condition edge cases passed
    pub boundary_condition_edge_cases: CheckStatus,
    /// Numerical stability edge cases passed
    pub numerical_stability_edge_cases: CheckStatus,
    /// Physical constraint validation passed
    pub physical_constraint_validation: CheckStatus,
    /// Convergence algorithm robustness passed
    pub convergence_algorithm_robustness: CheckStatus,
    /// Implementation edge cases passed
    pub implementation_edge_cases: CheckStatus,
}

/// Edge case testing results
#[derive(Debug, Clone)]
pub struct EdgeCaseTesting {
    /// Edge case validation results
    pub edge_case_results: Vec<EdgeCaseResult>,
}
