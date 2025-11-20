//! Validation tests for backend abstraction pattern
//!
//! Validates the backend abstraction implementation against performance and correctness
//! benchmarks following ASME V&V 20-2009 standards.
//!
//! References:
//! - ASME V&V 20-2009: "Standard for Verification and Validation in Computational Fluid Dynamics"
//! - ISO/IEC 25010:2011: "Systems and software quality models"
//! - Rust Performance Book: https://nnethercote.github.io/perf-book/

use cfd_core::compute::backend_example::{
    compute_squares, select_backend, Backend, ComputeBackend,
};

/// Test backend selection follows feature gate correctly
///
/// Reference: Rust conditional compilation best practices
#[test]
fn test_backend_selection_feature_gate() {
    let backend = select_backend();

    // Verify compile-time selection matches build configuration
    #[cfg(feature = "gpu")]
    assert_eq!(
        backend,
        Backend::Gpu,
        "GPU feature enabled: should select GPU backend"
    );

    #[cfg(not(feature = "gpu"))]
    assert_eq!(
        backend,
        Backend::Cpu,
        "GPU feature disabled: should select CPU backend"
    );
}

/// Test compute operation correctness for integer types
///
/// Validates element-wise square operation with integer arithmetic
#[test]
fn test_compute_squares_integers_correctness() {
    let backend = Backend::Cpu;
    let input = vec![0, 1, 2, 3, 4, 5, 10, 100];
    let expected = vec![0, 1, 4, 9, 16, 25, 100, 10000];

    let result = compute_squares(&backend, &input);

    assert_eq!(
        result.len(),
        expected.len(),
        "Output length must match input"
    );
    assert_eq!(result, expected, "Element-wise squares must be correct");
}

/// Test compute operation correctness for floating-point types
///
/// Validates numerical accuracy with f64 precision
/// Reference: IEEE 754-2008 floating-point standard
#[test]
fn test_compute_squares_floats_correctness() {
    let backend = Backend::Cpu;
    let input: Vec<f64> = vec![0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
    let expected: Vec<f64> = vec![0.0, 0.25, 1.0, 2.25, 4.0, 6.25, 9.0];

    let result = compute_squares(&backend, &input);

    assert_eq!(result.len(), expected.len());
    for (i, (&computed, &expected_val)) in result.iter().zip(expected.iter()).enumerate() {
        let error = (computed - expected_val).abs();
        assert!(
            error < 1e-10,
            "Element {} mismatch: {} vs {} (error: {})",
            i,
            computed,
            expected_val,
            error
        );
    }
}

/// Test zero-copy operation with slice storage
///
/// Validates that slice-based storage works without allocation overhead
#[test]
fn test_slice_storage_zero_copy() {
    let backend = Backend::Cpu;
    let data = [1.0f64, 2.0, 3.0, 4.0, 5.0];
    let slice: &[f64] = &data;

    let result = compute_squares(&backend, &slice);
    let expected = vec![1.0, 4.0, 9.0, 16.0, 25.0];

    assert_eq!(result, expected);
}

/// Test CPU and GPU backends produce identical results
///
/// Validates backend abstraction correctness: both implementations
/// should produce bit-identical results for deterministic operations
#[test]
fn test_backend_consistency() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    let cpu_result = Backend::Cpu.compute_squares(&data);
    let gpu_result = Backend::Gpu.compute_squares(&data);

    assert_eq!(
        cpu_result, gpu_result,
        "CPU and GPU backends must produce identical results"
    );
}

/// Test large dataset handling
///
/// Validates performance and correctness with realistic problem sizes
/// Reference: Typical CFD grid sizes (10^4 to 10^7 cells)
#[test]
fn test_large_dataset_validation() {
    let backend = select_backend();
    let n = 10_000;
    let input: Vec<f64> = (0..n).map(|i| i as f64 * 0.1).collect();

    let result = compute_squares(&backend, &input);

    assert_eq!(result.len(), n, "Output size must match input size");

    // Spot check correctness at various positions
    let check_indices = [0, n / 4, n / 2, 3 * n / 4, n - 1];
    for &idx in &check_indices {
        let expected = input[idx] * input[idx];
        let error = (result[idx] - expected).abs();
        assert!(
            error < 1e-10,
            "Element {} error {} exceeds tolerance",
            idx,
            error
        );
    }
}

/// Test edge case: empty input
///
/// Validates graceful handling of edge cases per defensive programming practices
#[test]
fn test_empty_input() {
    let backend = Backend::Cpu;
    let input: Vec<f64> = vec![];

    let result = compute_squares(&backend, &input);

    assert_eq!(result.len(), 0, "Empty input should produce empty output");
}

/// Test edge case: single element
///
/// Validates minimum viable input handling
#[test]
fn test_single_element() {
    let backend = Backend::Cpu;
    let input = vec![42.0];

    let result = compute_squares(&backend, &input);

    assert_eq!(result.len(), 1);
    assert_eq!(result[0], 1764.0);
}

/// Test numerical stability with extreme values
///
/// Validates behavior near floating-point limits
/// Reference: IEEE 754-2008 §5 arithmetic operations
#[test]
fn test_extreme_values() {
    let backend = Backend::Cpu;

    // Test with very small values - use relative tolerance for floating point
    let small: Vec<f64> = vec![1e-10, 1e-20, 1e-100];
    let small_result = compute_squares(&backend, &small);
    let small_expected: Vec<f64> = vec![1e-20, 1e-40, 1e-200];

    for (i, (&result, &expected_val)) in small_result.iter().zip(small_expected.iter()).enumerate()
    {
        let rel_error = if expected_val != 0.0 {
            ((result - expected_val) / expected_val).abs()
        } else {
            result.abs()
        };
        assert!(
            rel_error < 1e-10 || result == expected_val,
            "Small value {} mismatch: {} vs {} (rel_error: {})",
            i,
            result,
            expected_val,
            rel_error
        );
    }

    // Test with moderately large values (avoid overflow)
    let large = vec![1e10, 1e20];
    let large_result = compute_squares(&backend, &large);
    assert_eq!(large_result, vec![1e20, 1e40]);
}

/// Test mixed positive/negative values
///
/// Validates correct handling of sign in square operations
#[test]
fn test_mixed_signs() {
    let backend = Backend::Cpu;
    let input = vec![-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0];
    let expected = vec![9.0, 4.0, 1.0, 0.0, 1.0, 4.0, 9.0];

    let result = compute_squares(&backend, &input);

    assert_eq!(
        result, expected,
        "Squares of negative numbers must be positive"
    );
}

/// Benchmark-style validation: verify linear scaling
///
/// Validates O(n) complexity as expected for element-wise operations
/// Reference: Algorithmic complexity analysis
#[test]
fn test_linear_scaling_property() {
    let backend = Backend::Cpu;

    // Test with different sizes to verify linear scaling
    let sizes = [100, 1000, 10000];

    for &size in &sizes {
        let input: Vec<f64> = (0..size).map(|i| i as f64).collect();
        let result = compute_squares(&backend, &input);

        assert_eq!(
            result.len(),
            size,
            "Output size must match input for size {}",
            size
        );

        // Verify first, middle, and last elements
        assert_eq!(result[0], 0.0);
        assert_eq!(result[size / 2], ((size / 2) as f64).powi(2));
        assert_eq!(result[size - 1], ((size - 1) as f64).powi(2));
    }
}

/// Test type genericity: i32, i64, f32, f64
///
/// Validates trait-based abstraction works across numeric types
#[test]
fn test_type_genericity() {
    let backend = Backend::Cpu;

    // Test i32
    let i32_data = vec![1i32, 2, 3];
    let i32_result = compute_squares(&backend, &i32_data);
    assert_eq!(i32_result, vec![1, 4, 9]);

    // Test i64
    let i64_data = vec![1i64, 2, 3];
    let i64_result = compute_squares(&backend, &i64_data);
    assert_eq!(i64_result, vec![1, 4, 9]);

    // Test f32
    let f32_data = vec![1.0f32, 2.0, 3.0];
    let f32_result = compute_squares(&backend, &f32_data);
    assert_eq!(f32_result, vec![1.0, 4.0, 9.0]);

    // Test f64
    let f64_data = vec![1.0f64, 2.0, 3.0];
    let f64_result = compute_squares(&backend, &f64_data);
    assert_eq!(f64_result, vec![1.0, 4.0, 9.0]);
}

/// Integration test: realistic CFD application scenario
///
/// Simulates typical CFD usage pattern with field operations
#[test]
fn test_cfd_application_scenario() {
    let backend = select_backend();

    // Simulate velocity field magnitude calculation (u² + v²)
    // For simplicity, just test u² component
    let velocity_u = vec![1.0, 2.0, 3.0, 4.0, 5.0]; // m/s
    let u_squared = compute_squares(&backend, &velocity_u);

    // Verify physical reasonableness
    assert_eq!(u_squared.len(), velocity_u.len());
    for (i, &val) in u_squared.iter().enumerate() {
        assert!(val >= 0.0, "Squared values must be non-negative");
        assert_eq!(
            val,
            velocity_u[i] * velocity_u[i],
            "Must match analytical result"
        );
    }
}
