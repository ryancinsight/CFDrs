//! Benchmark validation tests for CFD benchmark problems
//!
//! Tests the validation logic for standard CFD benchmarks:
//! - Backward-facing step (Gartling 1990)
//! - Flow over cylinder (Schäfer & Turek 1996)

use cfd_validation::benchmarks::{
    BackwardFacingStep, Benchmark, BenchmarkConfig, BenchmarkResult, FlowOverCylinder,
};

#[test]
fn test_backward_facing_step_reference_solution() {
    // Create benchmark with standard expansion ratio (ER=2.0)
    let step = BackwardFacingStep::<f64>::new(
        1.0,  // step_height
        2.0,  // channel_height
        10.0, // channel_length
        1.0,  // inlet_velocity
    );

    // Reference solution should be available
    let reference = step.reference_solution();
    assert!(reference.is_some(), "Reference solution should exist");

    let reference = reference.unwrap();
    assert_eq!(reference.name, "Backward Facing Step (Reference)");
    assert_eq!(
        reference.values.len(),
        1,
        "Should have one value (reattachment length)"
    );

    // Reference reattachment length should be around 6.0 * step_height
    // based on Gartling (1990) and Armaly et al. (1983)
    let reattachment = reference.values[0];
    assert!(
        reattachment > 0.0 && reattachment < 20.0,
        "Reattachment length should be physically reasonable"
    );

    // For moderate Re (200-400), expect x_r/h ≈ 6.0
    let normalized_reattachment = reattachment / 1.0; // Divide by step_height
    assert!(
        (normalized_reattachment - 6.0).abs() < 0.1,
        "Normalized reattachment should be ~6.0 for moderate Re"
    );
}

#[test]
fn test_backward_facing_step_validation_success() {
    let step = BackwardFacingStep::<f64>::new(1.0, 2.0, 10.0, 1.0);

    // Create a result that should pass validation
    // Reattachment length within 30% of reference (6.0)
    let result = BenchmarkResult {
        name: "Test Result".to_string(),
        values: vec![5.5], // Within 30% of 6.0
        errors: vec![],
        convergence: vec![1e-5], // Converged
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = step.validate(&result).expect("Validation should not error");
    assert!(
        is_valid,
        "Result within tolerance should validate successfully"
    );
}

#[test]
fn test_backward_facing_step_validation_tolerance() {
    let step = BackwardFacingStep::<f64>::new(1.0, 2.0, 10.0, 1.0);

    // Test boundary of 30% tolerance
    // Reference is 6.0, so 30% tolerance means [4.2, 7.8]

    // Just inside lower bound
    let result_low = BenchmarkResult {
        name: "Low Bound".to_string(),
        values: vec![4.3], // Just above 4.2
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        step.validate(&result_low).unwrap(),
        "Value just inside lower tolerance should pass"
    );

    // Just inside upper bound
    let result_high = BenchmarkResult {
        name: "High Bound".to_string(),
        values: vec![7.7], // Just below 7.8
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        step.validate(&result_high).unwrap(),
        "Value just inside upper tolerance should pass"
    );

    // Just outside upper bound
    let result_too_high = BenchmarkResult {
        name: "Too High".to_string(),
        values: vec![8.0], // Outside tolerance
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        !step.validate(&result_too_high).unwrap(),
        "Value outside tolerance should fail"
    );
}

#[test]
fn test_backward_facing_step_validation_failure_no_convergence() {
    let step = BackwardFacingStep::<f64>::new(1.0, 2.0, 10.0, 1.0);

    // Result with poor convergence
    let result = BenchmarkResult {
        name: "No Convergence".to_string(),
        values: vec![6.0], // Good reattachment length
        errors: vec![],
        convergence: vec![0.1], // Did not converge
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = step.validate(&result).unwrap();
    assert!(
        !is_valid,
        "Result without convergence should fail validation"
    );
}

#[test]
fn test_backward_facing_step_validation_failure_unphysical() {
    let step = BackwardFacingStep::<f64>::new(1.0, 2.0, 10.0, 1.0);

    // Unphysically large reattachment length
    let result = BenchmarkResult {
        name: "Unphysical".to_string(),
        values: vec![25.0], // Way too large (> 20 * step_height)
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = step.validate(&result).unwrap();
    assert!(!is_valid, "Unphysical result should fail validation");
}

#[test]
fn test_flow_over_cylinder_reference_solution() {
    // Create benchmark with standard parameters (Re=20 steady flow)
    let cylinder = FlowOverCylinder::<f64>::new(
        0.1,         // diameter
        (2.2, 0.41), // domain (length, height)
        1.0,         // inlet_velocity
    );

    // Reference solution should be available
    let reference = cylinder.reference_solution();
    assert!(reference.is_some(), "Reference solution should exist");

    let reference = reference.unwrap();
    assert_eq!(
        reference.name,
        "Flow Over Cylinder (Schäfer & Turek 1996, Re=20)"
    );
    assert_eq!(reference.values.len(), 2, "Should have two values (Cd, Cl)");

    // Check reference values match Schäfer & Turek (1996)
    let cd = reference.values[0];
    let cl = reference.values[1];

    // Cd should be around 5.57 for Re=20
    assert!(
        (cd - 5.57).abs() < 0.01,
        "Reference Cd should match Schäfer & Turek (5.57 ± 0.01)"
    );

    // Cl should be near zero with small asymmetry (0.0106)
    assert!(
        cl.abs() < 0.02,
        "Reference Cl should be near zero for Re=20"
    );

    // Check that errors are provided
    assert_eq!(reference.errors.len(), 2, "Should have error estimates");
}

#[test]
fn test_flow_over_cylinder_validation_success() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Create a result that should pass validation
    // Cd within 5% of reference (5.57), Cl near zero
    let result = BenchmarkResult {
        name: "Test Result".to_string(),
        values: vec![5.50, 0.02], // Cd close to 5.57, Cl small
        errors: vec![],
        convergence: vec![1e-5], // Converged
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = cylinder
        .validate(&result)
        .expect("Validation should not error");
    assert!(
        is_valid,
        "Result within tolerance should validate successfully"
    );
}

#[test]
fn test_flow_over_cylinder_validation_cd_tolerance() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Test 5% tolerance on Cd
    // Reference is 5.57, so 5% tolerance means [5.29, 5.85]

    // Just inside lower bound
    let result_low = BenchmarkResult {
        name: "Low Cd".to_string(),
        values: vec![5.30, 0.01], // Just above 5.29
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        cylinder.validate(&result_low).unwrap(),
        "Cd just inside lower tolerance should pass"
    );

    // Just outside lower bound
    let result_too_low = BenchmarkResult {
        name: "Too Low Cd".to_string(),
        values: vec![5.20, 0.01], // Below 5.29
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        !cylinder.validate(&result_too_low).unwrap(),
        "Cd outside lower tolerance should fail"
    );
}

#[test]
fn test_flow_over_cylinder_validation_cl_tolerance() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Test absolute tolerance on Cl (should be < 0.1)
    let result_good_cl = BenchmarkResult {
        name: "Good Cl".to_string(),
        values: vec![5.57, 0.05], // Cl within tolerance
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        cylinder.validate(&result_good_cl).unwrap(),
        "Small Cl should pass validation"
    );

    // Cl too large (symmetry breaking)
    let result_bad_cl = BenchmarkResult {
        name: "Bad Cl".to_string(),
        values: vec![5.57, 0.15], // Cl > 0.1
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };
    assert!(
        !cylinder.validate(&result_bad_cl).unwrap(),
        "Large Cl should fail validation"
    );
}

#[test]
fn test_flow_over_cylinder_validation_failure_no_convergence() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Result with poor convergence
    let result = BenchmarkResult {
        name: "No Convergence".to_string(),
        values: vec![5.57, 0.01], // Good coefficients
        errors: vec![],
        convergence: vec![0.01], // Did not converge
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = cylinder.validate(&result).unwrap();
    assert!(
        !is_valid,
        "Result without convergence should fail validation"
    );
}

#[test]
fn test_flow_over_cylinder_validation_failure_unphysical_cd() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Unphysically large drag coefficient
    let result = BenchmarkResult {
        name: "Unphysical Cd".to_string(),
        values: vec![25.0, 0.01], // Way too large
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = cylinder.validate(&result).unwrap();
    assert!(
        !is_valid,
        "Unphysical drag coefficient should fail validation"
    );
}

#[test]
fn test_flow_over_cylinder_validation_failure_unphysical_cl() {
    let cylinder = FlowOverCylinder::<f64>::new(0.1, (2.2, 0.41), 1.0);

    // Unphysically large lift coefficient
    let result = BenchmarkResult {
        name: "Unphysical Cl".to_string(),
        values: vec![5.57, 8.0], // Way too large for steady flow
        errors: vec![],
        convergence: vec![1e-5],
        execution_time: 0.0,
        metadata: std::collections::HashMap::new(),
    };

    let is_valid = cylinder.validate(&result).unwrap();
    assert!(
        !is_valid,
        "Unphysical lift coefficient should fail validation"
    );
}

#[test]
fn test_benchmark_with_generic_types() {
    // Test that benchmarks work with f32 as well as f64
    let step_f32 = BackwardFacingStep::<f32>::new(1.0, 2.0, 10.0, 1.0);
    let reference_f32 = step_f32.reference_solution();
    assert!(reference_f32.is_some(), "Should work with f32");

    let cylinder_f32 = FlowOverCylinder::<f32>::new(0.1, (2.2, 0.41), 1.0);
    let reference_cylinder_f32 = cylinder_f32.reference_solution();
    assert!(reference_cylinder_f32.is_some(), "Should work with f32");
}

#[test]
fn test_benchmark_run_integration() {
    // Integration test: Run benchmark and validate results
    let step = BackwardFacingStep::<f64>::new(1.0, 2.0, 10.0, 1.0);

    let config = BenchmarkConfig::<f64>::default();

    // Run the benchmark
    let result = step.run(&config);
    assert!(result.is_ok(), "Benchmark should run without error");

    let result = result.unwrap();

    // Check that result has expected structure
    assert!(!result.values.is_empty(), "Should have reattachment length");
    assert!(
        !result.convergence.is_empty(),
        "Should have convergence history"
    );

    // Validate the result (may pass or fail depending on solver implementation)
    let is_valid = step.validate(&result);
    assert!(is_valid.is_ok(), "Validation should not error");
}
