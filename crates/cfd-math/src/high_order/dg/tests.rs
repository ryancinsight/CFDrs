use leto::{Array1, Array2};

/// The copy-then-accumulate pair this replaced, as the differential
/// oracle: allocate the column, then add it scaled.
fn add_assign_scaled_via_column_copy(
    target: &mut Array1<f64>,
    matrix: &Array2<f64>,
    col: usize,
    scale: f64,
) {
    let source = column_vector(matrix, col);
    for row in 0..target.shape()[0] {
        target[row] += scale * source[row];
    }
}

#[test]
fn column_accumulate_is_bitwise_identical_to_copying_the_column() {
    // Values chosen so the products are not exactly representable and a
    // reordered sum would show: a tolerance-based assertion here would
    // pass for a reordering that breaks the claim.
    let matrix = Array2::from_shape_fn([4, 3], |[r, c]| {
        0.1 * (r as f64 + 1.0) + 1e-17 * (c as f64 + 1.0)
    });
    let scales = [0.3, -7.125, 1e12, f64::MIN_POSITIVE];

    for (col, scale) in (0..3).zip(scales) {
        let mut in_place = Array1::from_shape_fn([4], |[r]| 1.0 / (r as f64 + 3.0));
        let mut via_copy = in_place.clone();
        vector_add_assign_scaled_column(&mut in_place, &matrix, col, scale);
        add_assign_scaled_via_column_copy(&mut via_copy, &matrix, col, scale);
        for row in 0..4 {
            assert_eq!(
                in_place[row].to_bits(),
                via_copy[row].to_bits(),
                "column {col} row {row} drifted from the copying form",
            );
        }
    }
}

#[test]
#[should_panic(expected = "requested column is in bounds")]
fn column_accumulate_rejects_a_column_past_the_matrix() {
    let matrix = Array2::zeros([2, 2]);
    let mut target = Array1::zeros([2]);
    vector_add_assign_scaled_column(&mut target, &matrix, 2, 1.0);
}

#[test]
#[should_panic(expected = "equal length")]
fn column_accumulate_rejects_a_target_of_the_wrong_length() {
    let matrix = Array2::zeros([3, 2]);
    let mut target = Array1::zeros([2]);
    vector_add_assign_scaled_column(&mut target, &matrix, 0, 1.0);
}

use super::*;
use eunomia::assert_relative_eq;

#[test]
fn test_dg_solution() {
    // Create a DG solution with a quadratic polynomial
    let order = 2;
    let num_components = 1;
    let mut sol = DGSolution::new(order, num_components).expect("expected value");

    // Set the coefficients for u(x) = 1 + x + x²
    // In the Legendre basis: 4/3*P₀(x) + P₁(x) + (2/3)*P₂(x)
    sol.coefficients[[0, 0]] = 4.0 / 3.0; // P₀ term
    sol.coefficients[[0, 1]] = 1.0; // P₁ term
    sol.coefficients[[0, 2]] = 2.0 / 3.0; // P₂ term

    // Test evaluation at x = 1.0
    let u1 = sol.evaluate(1.0);
    assert_relative_eq!(u1[0], 3.0, epsilon = 1e-10);

    // Test evaluation at x = 0.0
    let u0 = sol.evaluate(0.0);
    assert_relative_eq!(u0[0], 1.0, epsilon = 1e-10);

    // Test L² norm
    // ∫(1 + x + x²)² dx from -1 to 1 = 4.4
    let norm = sol.l2_norm();
    assert_relative_eq!(norm, (4.4f64).sqrt(), epsilon = 1e-10);
}

#[test]
fn test_dg_solution_average() {
    // Create a DG solution with a constant function
    let order = 2;
    let num_components = 1;
    let mut sol = DGSolution::new(order, num_components).expect("expected value");

    // Set the coefficients for u(x) = 2.0
    // Since P₀(x) = 1.0, u(x) = 2.0 * P₀(x)
    sol.coefficients[[0, 0]] = 2.0;

    // The average should be 2.0
    let avg = sol.average();
    assert_relative_eq!(avg[0], 2.0, epsilon = 1e-10);
}
