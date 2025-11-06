//! Automated validation suite combining MMS, Richardson extrapolation, and performance benchmarking
//!
//! This comprehensive test suite validates CFD implementations through:
//! - Method of Manufactured Solutions (MMS) for various physics
//! - Richardson extrapolation for convergence verification
//! - Performance benchmarking for efficiency assessment
//! - Automated regression detection and reporting

use cfd_validation::benchmarking::{BenchmarkSuite, BenchmarkConfig, ExportFormat};
use cfd_validation::manufactured::{ManufacturedDiffusion, ManufacturedNavierStokes, ManufacturedKEpsilon, TaylorGreenManufactured};
use cfd_validation::manufactured::richardson_integration::{MmsRichardsonStudy, RichardsonMmsResult};
use cfd_validation::manufactured::ManufacturedSolution;
use std::collections::HashMap;

/// Comprehensive MMS validation for different physics
#[test]
fn test_comprehensive_mms_validation() {
    println!("Running Comprehensive MMS Validation Suite...");

    // Test diffusion equation MMS
    let diffusion_mms = ManufacturedDiffusion::new(1.0, 1.0, 1.0);
    validate_mms_solution(&diffusion_mms, "Diffusion_2D");

    // Test Navier-Stokes MMS
    let ns_mms = ManufacturedNavierStokes::new(1.0, 1.0, 1.0, 0.01);
    validate_mms_solution(&ns_mms, "NavierStokes_2D");

    // Test Taylor-Green vortex MMS
    let tg_mms = TaylorGreenManufactured::new(1.0, 1.0, 1.0);
    validate_mms_solution(&tg_mms, "TaylorGreen_2D");

    // Test k-ε turbulence MMS
    let keps_mms = ManufacturedKEpsilon::<f64>::new(1.0, 1.0, 1.0, 0.01);
    validate_mms_solution(&keps_mms, "KEpsilon_Turbulence");

    println!("✓ Comprehensive MMS validation completed");
}

/// Richardson extrapolation validation for convergence verification
#[test]
fn test_richardson_extrapolation_validation() {
    println!("Running Richardson Extrapolation Validation...");

    // Test with different manufactured solutions
    let test_cases = vec![
        ("Diffusion", ManufacturedDiffusion::new(1.0, 1.0, 1.0)),
        ("TaylorGreen", TaylorGreenManufactured::new(1.0, 1.0, 1.0)),
        ("NavierStokes", ManufacturedNavierStokes::new(1.0, 1.0, 1.0, 0.01)),
    ];

    for (name, mms) in test_cases {
        println!("Testing Richardson extrapolation for {}...", name);

        // Create Richardson study with geometric refinement
        let study = MmsRichardsonStudy::with_geometric_refinement(
            Box::new(mms),
            Box::new(cfd_validation::geometry::RectangularDomain::new(0.0, 1.0, 0.0, 1.0)),
            4, // 4 grid levels for robust convergence analysis
            0.1, // base grid size
            1.0, // evaluation time
        ).unwrap();

        // Run the study
        let result = study.run_study().unwrap();

        // Validate convergence properties
        validate_convergence_study(&result, name);

        println!("  {}: Order={:.3}, R²={:.6}, Asymptotic={}",
                name,
                result.convergence_study.convergence_rate,
                result.convergence_study.r_squared,
                result.convergence_study.is_asymptotic());
    }

    println!("✓ Richardson extrapolation validation completed");
}

/// Performance benchmarking validation
#[test]
fn test_performance_benchmarking_validation() {
    println!("Running Performance Benchmarking Validation...");

    // Configure benchmark suite
    let config = BenchmarkConfig {
        iterations: 3, // Reduced for testing
        enable_memory: true,
        enable_scaling: true,
        detailed_reporting: false,
        ..Default::default()
    };

    let suite = BenchmarkSuite::with_config(config);

    // Run benchmark suite
    let results = suite.run_full_suite().unwrap();

    // Validate benchmark results
    assert!(!results.is_empty(), "Should have benchmark results");

    let passed_count = results.iter()
        .filter(|r| matches!(r.status, cfd_validation::benchmarking::suite::BenchmarkStatus::Passed))
        .count();

    let total_count = results.len();
    let pass_rate = passed_count as f64 / total_count as f64;

    println!("Benchmark Results: {}/{} passed ({:.1}%)",
            passed_count, total_count, pass_rate * 100.0);

    assert!(pass_rate > 0.8, "Benchmark pass rate should be >80%, got {:.1}%", pass_rate * 100.0);

    // Export results for analysis
    let csv_report = suite.export_results(&results, ExportFormat::Csv).unwrap();
    println!("CSV Report length: {} characters", csv_report.len());

    println!("✓ Performance benchmarking validation completed");
}

/// Automated regression detection validation
#[test]
fn test_regression_detection_validation() {
    println!("Running Regression Detection Validation...");

    // Create baseline data (simulated previous run)
    let mut baseline = HashMap::new();
    baseline.insert("CFD_Matrix_Operations".to_string(), 0.01); // 10ms baseline
    baseline.insert("CFD_Vector_Operations".to_string(), 0.005); // 5ms baseline

    let config = BenchmarkConfig {
        iterations: 2,
        enable_memory: false,
        enable_scaling: false,
        regression_threshold: 10.0, // 10% threshold
        baseline_data: Some(baseline),
        ..Default::default()
    };

    let suite = BenchmarkSuite::with_config(config);

    // Run performance benchmarks only
    let results = suite.run_performance_benchmarks().unwrap();

    // Check for regressions (this will depend on actual performance)
    let regression_count = results.iter()
        .filter(|r| r.regression_detected.is_some())
        .count();

    println!("Regression analysis: {}/{} benchmarks showed changes from baseline",
            regression_count, results.len());

    // Export regression report
    let regression_report = suite.export_results(&results, ExportFormat::Text).unwrap();
    println!("Regression Report:\n{}", regression_report);

    println!("✓ Regression detection validation completed");
}

/// Cross-validation between different verification methods
#[test]
fn test_cross_validation_methods() {
    println!("Running Cross-Validation Methods Test...");

    // Use the same manufactured solution with different validation approaches
    let mms = ManufacturedDiffusion::new(1.0, 1.0, 1.0);

    // Method 1: Direct MMS validation
    validate_mms_solution(&mms, "CrossVal_Direct");

    // Method 2: Richardson extrapolation
    let study = MmsRichardsonStudy::with_geometric_refinement(
        Box::new(mms),
        Box::new(cfd_validation::geometry::RectangularDomain::new(0.0, 1.0, 0.0, 1.0)),
        3,
        0.1,
        1.0,
    ).unwrap();

    let richardson_result = study.run_study().unwrap();

    // Method 3: Performance benchmarking (on MMS evaluation)
    let perf_result = benchmark_mms_evaluation(&mms);

    // Cross-validate results
    assert!(richardson_result.convergence_study.convergence_rate > 1.5,
           "Diffusion should show ~2nd order convergence, got {:.3}",
           richardson_result.convergence_study.convergence_rate);

    assert!(perf_result > 0.0, "MMS evaluation should take measurable time");

    println!("Cross-validation: Richardson order={:.3}, Performance={:.6}ms",
            richardson_result.convergence_study.convergence_rate, perf_result * 1000.0);

    println!("✓ Cross-validation methods test completed");
}

/// Helper function to validate MMS solution properties
fn validate_mms_solution<T: ManufacturedSolution<f64> + Clone>(mms: &T, test_name: &str) {
    // Test solution properties at multiple points
    let test_points = [
        (0.25, 0.25, 0.0, 0.5),
        (0.5, 0.5, 0.0, 1.0),
        (0.75, 0.75, 0.0, 1.5),
    ];

    for (x, y, z, t) in test_points {
        let solution = mms.exact_solution(x, y, z, t);
        let source = mms.source_term(x, y, z, t);

        // Basic sanity checks
        assert!(solution.is_finite(), "{}: Solution should be finite at ({}, {}, {}, {})",
               test_name, x, y, z, t);
        assert!(source.is_finite(), "{}: Source term should be finite at ({}, {}, {}, {})",
               test_name, x, y, z, t);

        // For manufactured solutions, solution should be bounded
        assert!(solution.abs() < 1000.0, "{}: Solution magnitude too large: {} at ({}, {}, {}, {})",
               test_name, solution, x, y, z, t);
    }

    println!("  {}: MMS validation passed", test_name);
}

/// Helper function to validate convergence study results
fn validate_convergence_study(result: &RichardsonMmsResult<f64>, test_name: &str) {
    // Check convergence rate is reasonable (between 1 and 3 for CFD methods)
    assert!((1.0..=3.0).contains(&result.convergence_study.convergence_rate),
           "{}: Convergence rate {:.3} outside reasonable range [1, 3]",
           test_name, result.convergence_study.convergence_rate);

    // Check R-squared indicates good fit
    assert!(result.convergence_study.r_squared > 0.9,
           "{}: Poor convergence fit R²={:.6}",
           test_name, result.convergence_study.r_squared);

    // Check Richardson extrapolation gives consistent results
    if let (Some(order1), Some(order2)) = (
        result.richardson_results.first().map(|(_, o)| *o),
        result.richardson_results.get(1).map(|(_, o)| *o)
    ) {
        let order_diff = (order1 - order2).abs();
        assert!(order_diff < 0.5, "{}: Richardson orders inconsistent: {:.3} vs {:.3}",
               test_name, order1, order2);
    }

    // Check GCI values are reasonable (should be small for converged solutions)
    for &gci in &result.gci_values {
        if gci > 0.0 {
            assert!(gci < 1.0, "{}: GCI value {:.3} too large (should be < 1.0)",
                   test_name, gci);
        }
    }
}

/// Helper function to benchmark MMS evaluation performance
fn benchmark_mms_evaluation(mms: &impl ManufacturedSolution<f64>) -> f64 {
    use std::time::Instant;

    let num_evaluations = 10000;
    let mut total_time = 0.0;

    // Warm-up
    for _ in 0..100 {
        let _ = mms.exact_solution(0.5, 0.5, 0.0, 1.0);
        let _ = mms.source_term(0.5, 0.5, 0.0, 1.0);
    }

    // Timed evaluation
    let start = Instant::now();
    for i in 0..num_evaluations {
        let x = (i % 100) as f64 * 0.01;
        let y = ((i / 100) % 100) as f64 * 0.01;
        let t = (i % 10) as f64 * 0.1;

        let _solution = mms.exact_solution(x, y, 0.0, t);
        let _source = mms.source_term(x, y, 0.0, t);
    }
    let elapsed = start.elapsed();

    elapsed.as_secs_f64() / num_evaluations as f64
}

/// Integration test combining all validation methods
#[test]
fn test_integrated_validation_pipeline() {
    println!("Running Integrated Validation Pipeline...");

    // Step 1: MMS validation
    let mms = ManufacturedDiffusion::new(2.0, 1.5, 1.0);
    validate_mms_solution(&mms, "Integrated_Diffusion");

    // Step 2: Richardson convergence study
    let study = MmsRichardsonStudy::with_geometric_refinement(
        Box::new(mms),
        Box::new(cfd_validation::geometry::RectangularDomain::new(0.0, 1.0, 0.0, 1.0)),
        3,
        0.05,
        1.0,
    ).unwrap();

    let convergence_result = study.run_study().unwrap();
    validate_convergence_study(&convergence_result, "Integrated_Convergence");

    // Step 3: Performance assessment
    let perf_time = benchmark_mms_evaluation(&ManufacturedDiffusion::new(2.0, 1.5, 1.0));

    // Step 4: Validation report
    println!("Integrated Validation Report:");
    println!("  MMS: ✓ Validated");
    println!("  Richardson: Order={:.3}, R²={:.6}",
            convergence_result.convergence_study.convergence_rate,
            convergence_result.convergence_study.r_squared);
    println!("  Performance: {:.3}μs per evaluation", perf_time * 1_000_000.0);
    println!("  Overall: ✓ PASSED");

    assert!(convergence_result.convergence_study.is_asymptotic(),
           "Solution should be in asymptotic convergence range");

    println!("✓ Integrated validation pipeline completed");
}

/// Property-based tests for validation robustness
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Test MMS validation with various parameters
        #[test]
        fn test_mms_parameter_robustness(
            kx in 0.1f64..5.0,
            ky in 0.1f64..5.0,
            amplitude in 0.1f64..2.0,
            nu in 1e-5f64..0.1
        ) {
            let mms = ManufacturedDiffusion::new(kx, ky, amplitude);
            let x = 0.5;
            let y = 0.5;
            let t = 1.0;

            let solution = mms.exact_solution(x, y, 0.0, t);
            let source = mms.source_term(x, y, 0.0, t);

            prop_assert!(solution.is_finite(), "Solution must be finite");
            prop_assert!(source.is_finite(), "Source must be finite");
            prop_assert!(solution.abs() < 1000.0, "Solution magnitude reasonable");
        }

        /// Test Richardson study with various configurations
        #[test]
        fn test_richardson_study_robustness(
            base_size in 0.01f64..0.5,
            eval_time in 0.1f64..2.0,
            levels in 2usize..5
        ) {
            let mms = ManufacturedDiffusion::new(1.0, 1.0, 1.0);
            let geometry = cfd_validation::geometry::RectangularDomain::new(0.0, 1.0, 0.0, 1.0);

            let study = MmsRichardsonStudy::with_geometric_refinement(
                Box::new(mms),
                Box::new(geometry),
                levels,
                base_size,
                eval_time,
            );

            prop_assert!(study.is_ok(), "Richardson study should be constructible");
            let study = study.unwrap();

            let result = study.run_study();
            prop_assert!(result.is_ok(), "Richardson study should run successfully");
        }
    }
}

