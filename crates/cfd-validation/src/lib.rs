//! Validation framework for CFD simulations
//!
//! This module provides tools for validating CFD implementations against:
//! - Analytical solutions (Poiseuille, Couette, Womersley)
//! - Benchmark problems (Lid-driven cavity, Flow over cylinder)
//! - Literature results (Sangalli et al., Bharadvaj et al.)
//! - Method of manufactured solutions (MMS)
//!
//! ## Mathematical Foundations
//!
//! ### 1. Non-Newtonian Rheology: Casson Model
//! The Casson model is widely used for blood flow, accounting for yield stress $(\tau_y)$ and shear-thinning behavior:
//! $$\sqrt{\tau} = \sqrt{\tau_y} + \sqrt{\mu_p \dot{\gamma}}$$
//! for $\tau > \tau_y$, where $\mu_p$ is the plastic viscosity and $\dot{\gamma}$ is the shear rate.
//!
//! ### 2. Murray's Law for Vascular Branching
//! For optimal metabolic efficiency in vascular networks, the parent diameter $D_0$ and daughter diameters $D_i$ must satisfy:
//! $D_0^3 = \sum D_i^3$
//! This library validates branching geometries against this cubic relationship.
//!
//! ### 3. Richardson Extrapolation and GCI
//! To verify grid convergence, we compute the Grid Convergence Index (GCI):
//! $GCI = \frac{1.25 |\epsilon|}{r^p - 1}$
//! where $\epsilon$ is the relative error between grids, $r$ is the refinement ratio, and $p$ is the observed order of accuracy.
//!
//! ### 4. Bernoulli and Cavitation Number
//! Venturi flow validation utilizes the Bernoulli principle and the cavitation number $\sigma$:
//! $\sigma = \frac{p_\infty - p_v}{\frac{1}{2}\rho v_\infty^2}$
//! determining the inception of vapor phase transition at the throat.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// CFD validation allows
#![allow(clippy::similar_names)] // Test variables often similar (u1,u2; err1,err2)
#![allow(clippy::cast_precision_loss)] // Acceptable in validation calculations
#![allow(clippy::must_use_candidate)] // Validation utilities often used in expressions
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::ptr_arg)] // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations
#![allow(clippy::approx_constant)] // Fallback constants for generic numerical types
#![allow(clippy::too_many_lines)] // Complex validation/benchmark functions need detailed implementation
#![allow(clippy::needless_range_loop)] // Explicit indexing clearer for multi-dimensional CFD arrays
#![allow(clippy::used_underscore_binding)] // Underscore prefixed bindings used for intentional partial use

pub mod adaptive_mesh;
pub mod algorithm_complexity;
pub mod analytical;
pub mod analytical_benchmarks;
pub mod benchmarking;
pub mod benchmarks;
pub mod conservation;
pub mod convergence;
pub mod edge_case_testing;
pub mod error_metrics;
pub mod geometry;
pub mod literature;
pub mod manufactured;
pub mod numerical;
pub mod reporting;
pub mod solutions;
pub mod time_integration;

// The public modules are the primary API.
// Users should access types through the module hierarchy:
//   use cfd_validation::analytical::AnalyticalSolution;
//   use cfd_validation::convergence::GridConvergenceIndex;
//   use cfd_validation::conservation::ConservationChecker;
//
// This hierarchical structure provides clear organization:
// - analytical: Known exact solutions (Poiseuille, Couette, etc.)
// - benchmarks: Standard test cases (lid-driven cavity, etc.)
// - conservation: Mass/momentum/energy conservation checks
// - convergence: Grid/temporal convergence studies
// - error_metrics: L2 norm, L∞ norm, etc.
// - literature: Comparison with published results
// - manufactured: Method of manufactured solutions
// - numerical: Numerical analysis tools
// - solutions: Solution comparison utilities
// - time_integration: Temporal accuracy analysis

/// Run comprehensive performance profiling for CFD algorithms
///
/// This function performs detailed performance analysis including:
/// - Algorithm complexity documentation (Big-O analysis)
/// - Memory bandwidth and cache efficiency measurements
/// - Performance recommendations based on empirical data
/// - Literature-backed complexity analysis
///
/// # Returns
/// Performance profiling report with detailed metrics and recommendations
///
/// # Examples
/// ```rust,no_run
/// use cfd_validation::run_performance_profiling;
///
/// let report = run_performance_profiling().expect("Performance profiling failed");
/// println!("Total performance: {:.2} GFLOPS", report.summary.total_gflops);
/// ```
pub fn run_performance_profiling(
) -> cfd_core::error::Result<benchmarking::production::PerformanceReport> {
    use benchmarking::production::PerformanceProfiler;

    let profiler = PerformanceProfiler::new();
    profiler.run_comprehensive_profiling()
}

/// Run comprehensive stability analysis for CFD time-stepping schemes
///
/// This function performs detailed stability analysis including:
/// - Stability region computation for Runge-Kutta methods
/// - CFL condition verification for various flow regimes
/// - Von Neumann stability analysis for linear PDEs
/// - Stability recommendations based on literature standards
///
/// # Returns
/// Stability analysis report with detailed stability metrics and recommendations
///
/// # Examples
/// ```rust,no_run
/// use cfd_validation::run_stability_analysis;
///
/// let report = run_stability_analysis().expect("Stability analysis failed");
/// println!("Overall stability score: {:.1}%", report.overall_assessment.overall_score * 100.0);
/// ```
pub fn run_stability_analysis(
) -> cfd_core::error::Result<time_integration::stability_analysis::StabilityAnalysisReport<f64>> {
    use time_integration::stability_analysis::StabilityAnalysisRunner;

    let analyzer = StabilityAnalysisRunner::<f64>::new();
    analyzer.run_comprehensive_stability_analysis()
}

/// Get the comprehensive CFD algorithm complexity registry
///
/// This function provides access to detailed Big-O complexity analysis for all
/// major CFD algorithms implemented in the codebase, including time/space complexity,
/// memory access patterns, cache efficiency, and parallel scalability metrics.
///
/// # Returns
/// Algorithm complexity registry with detailed performance analysis
///
/// # Examples
/// ```rust
/// use cfd_validation::get_algorithm_complexity_registry;
///
/// let registry = get_algorithm_complexity_registry();
/// if let Some(cg_info) = registry.get("ConjugateGradient") {
///     println!("CG Time Complexity: {}", cg_info.time_complexity);
///     println!("Cache Efficiency: {:.1}%", cg_info.cache_efficiency * 100.0);
/// }
/// ```
pub fn get_algorithm_complexity_registry() -> algorithm_complexity::AlgorithmComplexityRegistry {
    algorithm_complexity::AlgorithmComplexityRegistry::new()
}

/// Run comprehensive edge case testing for CFD algorithms
///
/// This function performs thorough validation of boundary conditions, numerical stability,
/// convergence algorithms, physical constraints, and implementation edge cases that
/// can cause failures in CFD simulations.
///
/// # Returns
/// Edge case test report with detailed validation results and recommendations
///
/// # Examples
/// ```rust,no_run
/// use cfd_validation::run_edge_case_testing;
///
/// let report = run_edge_case_testing().expect("Edge case testing failed");
/// println!("Overall pass rate: {:.1}%", report.overall_assessment.pass_rate * 100.0);
/// if report.overall_assessment.critical_failures > 0 {
///     println!("Warning: {} critical failures detected!", report.overall_assessment.critical_failures);
/// }
/// ```
pub fn run_edge_case_testing() -> cfd_core::error::Result<edge_case_testing::EdgeCaseReport> {
    use edge_case_testing::EdgeCaseTestSuite;

    let test_suite: EdgeCaseTestSuite<f64> = EdgeCaseTestSuite::new();
    test_suite.run_comprehensive_edge_case_tests()
}
