use super::super::DGSolution;

use super::*;
use crate::high_order::dg::matrix_from_vec;
use eunomia::assert_relative_eq;

#[test]
fn test_minmod_limiter() {
    let mut solution = DGSolution::new(2, 1).expect("expected value");
    solution.coefficients = matrix_from_vec(1, 3, vec![1.0, 1.0, 0.5]);

    let mut left = DGSolution::new(2, 1).expect("expected value");
    left.coefficients = matrix_from_vec(1, 3, vec![0.0, 0.0, 0.0]);

    let mut right = DGSolution::new(2, 1).expect("expected value");
    right.coefficients = matrix_from_vec(1, 3, vec![2.0, 0.0, 0.0]);

    let limiter = MinmodLimiter;
    let mut params = LimiterParams::new(LimiterType::Minmod);
    params.adaptive = false;

    limiter
        .limit(&mut solution, &[left, right], &params)
        .expect("expected value");

    assert_eq!(solution.coefficients[[0, 0]], 1.0);
    assert_eq!(solution.coefficients[[0, 1]], 1.0);
    assert_eq!(solution.coefficients[[0, 2]], 0.0);
}

#[test]
fn test_tvb_limiter() {
    let mut solution = DGSolution::new(2, 1).expect("expected value");
    solution.coefficients = matrix_from_vec(1, 3, vec![1.0, 0.1, 0.01]);

    let mut left = DGSolution::new(2, 1).expect("expected value");
    left.coefficients = matrix_from_vec(1, 3, vec![0.9, 0.1, 0.0]);

    let mut right = DGSolution::new(2, 1).expect("expected value");
    right.coefficients = matrix_from_vec(1, 3, vec![1.1, 0.1, 0.0]);

    let limiter = TVBLimiter;
    let params = LimiterParams::new(LimiterType::TVB).with_tvb_m(1.0);

    limiter
        .limit(&mut solution, &[left, right], &params)
        .expect("expected value");

    assert_relative_eq!(solution.coefficients[[0, 0]], 1.0, epsilon = 1e-10);
    assert_relative_eq!(solution.coefficients[[0, 1]], 0.1, epsilon = 1e-10);
    assert_relative_eq!(solution.coefficients[[0, 2]], 0.01, epsilon = 1e-10);
}

#[test]
fn test_moment_limiter() {
    let mut solution = DGSolution::new(3, 1).expect("expected value");
    solution.coefficients = matrix_from_vec(1, 4, vec![1.0, 1.0, 0.5, 0.1]);

    let mut left = DGSolution::new(3, 1).expect("expected value");
    left.coefficients = matrix_from_vec(1, 4, vec![0.0, 0.0, 0.0, 0.0]);

    let mut right = DGSolution::new(3, 1).expect("expected value");
    right.coefficients = matrix_from_vec(1, 4, vec![2.0, 0.0, 0.0, 0.0]);

    let limiter = MomentLimiter;
    let mut params = LimiterParams::new(LimiterType::Moment);
    params.adaptive = false;

    limiter
        .limit(&mut solution, &[left, right], &params)
        .expect("expected value");

    assert_eq!(solution.coefficients[[0, 0]], 1.0);
    assert_eq!(solution.coefficients[[0, 1]], 1.0);
    assert_eq!(solution.coefficients[[0, 2]], 0.0);
    assert_eq!(solution.coefficients[[0, 3]], 0.0);
}
