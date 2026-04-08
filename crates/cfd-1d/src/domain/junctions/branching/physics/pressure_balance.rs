//! Shared scalar pressure-balance utilities for branching junction solvers.

use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;

/// Scalar tolerances used by bracketed bisection solves.
#[derive(Debug, Clone, Copy)]
pub(super) struct ScalarSolveTolerances<T: RealField + Copy> {
    /// Absolute/relative target tolerance.
    pub value_tolerance: T,
    /// Interval-width tolerance.
    pub interval_tolerance: T,
    /// Maximum bisection iterations.
    pub max_iterations: usize,
}

impl<T: RealField + Copy + SafeFromF64> ScalarSolveTolerances<T> {
    /// Tolerances for flow-root solves over `[0, q_parent]`.
    pub(super) fn for_flow_interval(q_parent: T) -> Self {
        let tiny = T::from_f64_or_one(1e-18);
        Self {
            value_tolerance: T::from_f64_or_one(1e-12),
            interval_tolerance: q_parent.abs() * T::from_f64_or_one(1e-10) + tiny,
            max_iterations: 80,
        }
    }

    /// Tolerances for monotone target inversion.
    pub(super) fn for_target(target: T, interval_scale: T) -> Self {
        let tiny = T::from_f64_or_one(1e-18);
        Self {
            value_tolerance: target.abs() * T::from_f64_or_one(1e-10) + tiny,
            interval_tolerance: interval_scale.abs() * T::from_f64_or_one(1e-10) + tiny,
            max_iterations: 80,
        }
    }
}

/// Solve a bracketed scalar root with optional seed narrowing.
pub(super) fn bisect_root<T, F>(
    lower: T,
    upper: T,
    seed: Option<T>,
    tolerances: ScalarSolveTolerances<T>,
    residual: F,
) -> T
where
    T: RealField + Copy + SafeFromF64,
    F: Fn(T) -> T,
{
    let mut lower = lower;
    let mut upper = upper;

    if let Some(seed) = seed {
        let seed_residual = residual(seed);
        if seed_residual.abs() <= tolerances.value_tolerance {
            return seed;
        }
        if seed_residual < T::zero() {
            lower = seed;
        } else {
            upper = seed;
        }
    }

    for _ in 0..tolerances.max_iterations {
        let midpoint = (lower + upper) / T::from_f64_or_one(2.0);
        let midpoint_residual = residual(midpoint);

        if midpoint_residual.abs() <= tolerances.value_tolerance
            || (upper - lower).abs() <= tolerances.interval_tolerance
        {
            return midpoint;
        }

        if midpoint_residual < T::zero() {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }

    (lower + upper) / T::from_f64_or_one(2.0)
}

/// Invert a monotone nondecreasing scalar function by bracketed bisection.
pub(super) fn bisect_monotone_target<T, F>(
    lower: T,
    upper: T,
    target: T,
    tolerances: ScalarSolveTolerances<T>,
    value: F,
) -> T
where
    T: RealField + Copy + SafeFromF64,
    F: Fn(T) -> T,
{
    let mut lower = lower;
    let mut upper = upper;

    if value(upper) <= target {
        return upper;
    }

    for _ in 0..tolerances.max_iterations {
        let midpoint = (lower + upper) / T::from_f64_or_one(2.0);
        let midpoint_value = value(midpoint);

        if (midpoint_value - target).abs() <= tolerances.value_tolerance
            || (upper - lower).abs() <= tolerances.interval_tolerance
        {
            return midpoint;
        }

        if midpoint_value < target {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }

    (lower + upper) / T::from_f64_or_one(2.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bisect_root_finds_linear_zero() {
        let tolerances = ScalarSolveTolerances::for_flow_interval(1.0_f64);
        let root = bisect_root(0.0, 1.0, Some(0.8), tolerances, |x| x - 0.25);
        assert_relative_eq!(root, 0.25, epsilon = 1e-12);
    }

    #[test]
    fn test_bisect_monotone_target_inverts_linear_map() {
        let tolerances = ScalarSolveTolerances::for_target(0.75_f64, 1.0_f64);
        let value = bisect_monotone_target(0.0, 1.0, 0.75, tolerances, |x| x);
        assert_relative_eq!(value, 0.75, epsilon = 1e-12);
    }
}