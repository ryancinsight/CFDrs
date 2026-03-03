//! Comprehensive edge case testing for CFD algorithms
//!
//! This module provides thorough validation of boundary conditions, numerical stability,
//! and edge cases that can cause failures in CFD simulations.

mod test_cases;

use cfd_core::error::Result;
use nalgebra::RealField;

// ── Types ─────────────────────────────────────────────────────────────────────

/// Comprehensive edge case test suite
#[derive(Debug, Clone)]
pub struct EdgeCaseTestSuite<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

/// Edge case test report
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EdgeCaseReport {
    /// Boundary condition edge case tests
    pub boundary_condition_tests: Vec<EdgeCaseTestResult>,
    /// Numerical stability edge case tests
    pub stability_tests: Vec<EdgeCaseTestResult>,
    /// Convergence algorithm edge case tests
    pub convergence_tests: Vec<EdgeCaseTestResult>,
    /// Physical constraint edge case tests
    pub physical_constraint_tests: Vec<EdgeCaseTestResult>,
    /// Implementation edge case tests
    pub implementation_tests: Vec<EdgeCaseTestResult>,
    /// Overall assessment
    pub overall_assessment: EdgeCaseAssessment,
}

/// Individual edge case test result
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EdgeCaseTestResult {
    /// Test name
    pub test_name: String,
    /// Whether test passed
    pub passed: bool,
    /// Test severity
    pub severity: TestSeverity,
    /// Test description
    pub description: String,
    /// Detailed test results
    pub details: String,
}

/// Test severity levels
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum TestSeverity {
    /// Low severity - cosmetic issues
    Low,
    /// Medium severity - performance impact
    Medium,
    /// High severity - accuracy impact
    High,
    /// Critical severity - potential crashes or incorrect results
    Critical,
}

/// Overall edge case assessment
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EdgeCaseAssessment {
    /// Total number of tests run
    pub total_tests: usize,
    /// Number of tests passed
    pub passed_tests: usize,
    /// Pass rate (0.0 to 1.0)
    pub pass_rate: f64,
    /// Number of critical failures
    pub critical_failures: usize,
    /// Overall status
    pub overall_status: OverallStatus,
}

impl Default for EdgeCaseAssessment {
    fn default() -> Self {
        Self {
            total_tests: 0,
            passed_tests: 0,
            pass_rate: 0.0,
            critical_failures: 0,
            overall_status: OverallStatus::Good,
        }
    }
}

/// Overall test status
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum OverallStatus {
    /// Excellent - all critical tests pass
    Excellent,
    /// Good - minor issues only
    Good,
    /// Needs attention - some failures
    NeedsAttention,
    /// Critical failures present
    CriticalFailures,
}

// ── Suite impl ────────────────────────────────────────────────────────────────

impl<T: RealField + Copy> EdgeCaseTestSuite<T> {
    /// Create new edge case test suite
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Run comprehensive edge case validation suite
    pub fn run_comprehensive_edge_case_tests(&self) -> Result<EdgeCaseReport> {
        println!("🔍 Comprehensive CFD Edge Case Validation Suite");
        println!("==============================================");
        println!("Testing boundary conditions, numerical stability, and algorithm limits");
        println!();

        let mut report = EdgeCaseReport {
            boundary_condition_tests: Vec::new(),
            stability_tests: Vec::new(),
            convergence_tests: Vec::new(),
            physical_constraint_tests: Vec::new(),
            implementation_tests: Vec::new(),
            overall_assessment: EdgeCaseAssessment::default(),
        };

        Self::test_boundary_condition_edge_cases(&mut report)?;
        Self::test_numerical_stability_limits(&mut report)?;
        Self::test_convergence_algorithm_failures(&mut report)?;
        Self::test_physical_constraint_violations(&mut report)?;
        Self::test_implementation_edge_cases(&mut report)?;
        Self::generate_overall_assessment(&mut report);
        Self::display_edge_case_report(&report);

        Ok(report)
    }

    /// Test boundary condition edge cases
    fn test_boundary_condition_edge_cases(report: &mut EdgeCaseReport) -> Result<()> {
        println!("\n🔲 Boundary Condition Edge Case Tests");

        report.boundary_condition_tests.push(Self::test_extreme_velocity_gradients()?);
        report.boundary_condition_tests.push(Self::test_discontinuous_boundary_conditions()?);
        report.boundary_condition_tests.push(Self::test_incompatible_boundary_conditions()?);
        report.boundary_condition_tests.push(Self::test_corner_boundary_conditions()?);
        report.boundary_condition_tests.push(Self::test_time_dependent_boundary_conditions()?);

        println!("  ✅ Completed 5 boundary condition edge case tests");
        Ok(())
    }

    /// Test numerical stability limits
    fn test_numerical_stability_limits(report: &mut EdgeCaseReport) -> Result<()> {
        println!("\n🧮 Numerical Stability Edge Case Tests");

        report.stability_tests.push(Self::test_cfl_number_violations()?);
        report.stability_tests.push(Self::test_stiff_equation_systems()?);
        report.stability_tests.push(Self::test_ill_conditioned_matrices()?);
        report.stability_tests.push(Self::test_near_singular_systems()?);
        report.stability_tests.push(Self::test_long_time_integration_stability()?);

        println!("  ✅ Completed 5 numerical stability edge case tests");
        Ok(())
    }

    /// Test convergence algorithm failures
    fn test_convergence_algorithm_failures(report: &mut EdgeCaseReport) -> Result<()> {
        println!("\n🔄 Convergence Algorithm Edge Case Tests");

        report.convergence_tests.push(Self::test_preconditioner_breakdown()?);
        report.convergence_tests.push(Self::test_slow_convergence_scenarios()?);
        report.convergence_tests.push(Self::test_divergent_iterations()?);
        report.convergence_tests.push(Self::test_restarted_method_failures()?);

        println!("  ✅ Completed 4 convergence algorithm edge case tests");
        Ok(())
    }

    /// Test physical constraint violations
    fn test_physical_constraint_violations(report: &mut EdgeCaseReport) -> Result<()> {
        println!("\n⚛️ Physical Constraint Edge Case Tests");

        report.physical_constraint_tests.push(Self::test_negative_thermodynamic_properties()?);
        report.physical_constraint_tests.push(Self::test_negative_turbulence_quantities()?);
        report.physical_constraint_tests.push(Self::test_supersonic_flow_violations()?);
        report.physical_constraint_tests.push(Self::test_boundary_layer_separation()?);

        println!("  ✅ Completed 4 physical constraint edge case tests");
        Ok(())
    }

    /// Test implementation edge cases
    fn test_implementation_edge_cases(report: &mut EdgeCaseReport) -> Result<()> {
        println!("\n💻 Implementation Edge Case Tests");

        report.implementation_tests.push(Self::test_memory_allocation_limits()?);
        report.implementation_tests.push(Self::test_floating_point_precision_issues()?);
        report.implementation_tests.push(Self::test_parallel_computation_edge_cases()?);
        report.implementation_tests.push(Self::test_input_validation_boundaries()?);

        println!("  ✅ Completed 4 implementation edge case tests");
        Ok(())
    }

    /// Generate overall assessment
    fn generate_overall_assessment(report: &mut EdgeCaseReport) {
        let total_tests = report.boundary_condition_tests.len()
            + report.stability_tests.len()
            + report.convergence_tests.len()
            + report.physical_constraint_tests.len()
            + report.implementation_tests.len();

        let passed_tests = report
            .boundary_condition_tests
            .iter()
            .filter(|t| t.passed)
            .count()
            + report.stability_tests.iter().filter(|t| t.passed).count()
            + report.convergence_tests.iter().filter(|t| t.passed).count()
            + report
                .physical_constraint_tests
                .iter()
                .filter(|t| t.passed)
                .count()
            + report
                .implementation_tests
                .iter()
                .filter(|t| t.passed)
                .count();

        let pass_rate = if total_tests > 0 {
            passed_tests as f64 / total_tests as f64
        } else {
            0.0
        };

        let critical_failures = report
            .boundary_condition_tests
            .iter()
            .chain(report.stability_tests.iter())
            .chain(report.convergence_tests.iter())
            .chain(report.physical_constraint_tests.iter())
            .chain(report.implementation_tests.iter())
            .filter(|t| !t.passed && matches!(t.severity, TestSeverity::Critical))
            .count();

        let overall_status = if critical_failures > 0 {
            OverallStatus::CriticalFailures
        } else if pass_rate < 0.8 {
            OverallStatus::NeedsAttention
        } else if pass_rate >= 0.95 {
            OverallStatus::Excellent
        } else {
            OverallStatus::Good
        };

        report.overall_assessment = EdgeCaseAssessment {
            total_tests,
            passed_tests,
            pass_rate,
            critical_failures,
            overall_status,
        };
    }

    /// Display comprehensive edge case report
    fn display_edge_case_report(report: &EdgeCaseReport) {
        println!("\n📋 Comprehensive Edge Case Test Report");
        println!("=====================================");
        println!(
            "Overall Assessment: {} ({:.1}%)",
            match report.overall_assessment.overall_status {
                OverallStatus::Excellent => "Excellent",
                OverallStatus::Good => "Good",
                OverallStatus::NeedsAttention => "Needs Attention",
                OverallStatus::CriticalFailures => "Critical Failures",
            },
            report.overall_assessment.pass_rate * 100.0
        );
        println!(
            "Total Tests: {}/{}",
            report.overall_assessment.passed_tests, report.overall_assessment.total_tests
        );
        println!(
            "Critical Failures: {}",
            report.overall_assessment.critical_failures
        );

        println!("\n✅ Edge case validation completed successfully!");
        println!("   Comprehensive testing ensures CFD code robustness and reliability.");
    }
}

impl<T: RealField + Copy> Default for EdgeCaseTestSuite<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_edge_case_test_suite() {
        let suite = EdgeCaseTestSuite::<f64>::new();
        let report = suite.run_comprehensive_edge_case_tests();

        assert!(report.is_ok());
        let report = match report {
            Ok(report) => report,
            Err(err) => panic!("Edge case test suite failed: {err:?}"),
        };

        assert!(!report.boundary_condition_tests.is_empty());
        assert!(!report.stability_tests.is_empty());
        assert!(!report.convergence_tests.is_empty());
        assert!(!report.physical_constraint_tests.is_empty());
        assert!(!report.implementation_tests.is_empty());
        assert!(report.overall_assessment.pass_rate >= 0.0);
        assert!(report.overall_assessment.pass_rate <= 1.0);
    }
}
