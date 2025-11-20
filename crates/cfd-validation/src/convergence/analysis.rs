//! Convergence analysis utilities
//!
//! Provides tools for analyzing and classifying convergence behavior.

use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Convergence order classification
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceOrder<T: RealField + Copy> {
    /// Sub-linear convergence (p < 1)
    SubLinear,
    /// First-order convergence (p ≈ 1)
    FirstOrder,
    /// Super-linear convergence (1 < p < 2)
    SuperLinear,
    /// Second-order convergence (p ≈ 2)
    SecondOrder,
    /// Third-order convergence (p ≈ 3)
    ThirdOrder,
    /// Fourth-order convergence (p ≈ 4)
    FourthOrder,
    /// Spectral/exponential convergence
    Spectral,
    /// Custom order
    Custom(T),
}

impl<T: RealField + Copy + FromPrimitive> ConvergenceOrder<T> {
    /// Determine convergence order from observed rate
    pub fn from_rate(rate: T) -> Self {
        let tolerance = T::from_f64_or_one(0.15);

        if rate < T::from_f64(0.5).unwrap_or_else(T::zero) {
            Self::SubLinear
        } else if (rate - T::one()).abs() < tolerance {
            Self::FirstOrder
        } else if (rate - T::from_f64_or_one(2.0)).abs() < tolerance {
            Self::SecondOrder
        } else if (rate - T::from_f64_or_one(3.0)).abs() < tolerance {
            Self::ThirdOrder
        } else if (rate - T::from_f64_or_one(4.0)).abs() < tolerance {
            Self::FourthOrder
        } else if rate > T::from_f64_or_one(6.0) {
            Self::Spectral
        } else if rate > T::one() && rate < T::from_f64_or_one(2.0) {
            Self::SuperLinear
        } else {
            Self::Custom(rate)
        }
    }

    /// Get expected rate for this order
    pub fn expected_rate(&self) -> T {
        match self {
            Self::SubLinear => T::from_f64_or_one(0.5),
            Self::FirstOrder => T::one(),
            Self::SuperLinear => T::from_f64_or_one(1.5),
            Self::SecondOrder => T::from_f64_or_one(2.0),
            Self::ThirdOrder => T::from_f64_or_one(3.0),
            Self::FourthOrder => T::from_f64_or_one(4.0),
            Self::Spectral => T::from_f64_or_one(10.0),
            Self::Custom(rate) => *rate,
        }
    }

    /// Check if observed rate matches expected within tolerance
    pub fn matches(&self, observed_rate: T) -> bool {
        let expected = self.expected_rate();
        let tolerance = T::from_f64_or_one(0.25); // 25% tolerance
        (observed_rate - expected).abs() < expected * tolerance
    }
}

/// Convergence analysis tools
pub struct ConvergenceAnalysis;

impl ConvergenceAnalysis {
    /// Compute refinement ratio between consecutive grids
    pub fn refinement_ratio<T>(grid_sizes: &[T]) -> Vec<T>
    where
        T: RealField + Copy,
    {
        grid_sizes.windows(2).map(|w| w[0] / w[1]).collect()
    }

    /// Check if refinement is uniform (constant ratio)
    pub fn is_uniform_refinement<T>(grid_sizes: &[T]) -> bool
    where
        T: RealField + Copy + FromPrimitive + std::iter::Sum,
    {
        let ratios = Self::refinement_ratio(grid_sizes);
        if ratios.is_empty() {
            return true;
        }

        let mean_ratio =
            ratios.iter().copied().sum::<T>() / T::from_usize(ratios.len()).unwrap_or_else(T::one);
        let tolerance = T::from_f64_or_one(0.05); // 5% variation allowed

        ratios
            .iter()
            .all(|r| ((*r - mean_ratio) / mean_ratio).abs() < tolerance)
    }

    /// Estimate required grid size for target accuracy
    pub fn estimate_grid_for_accuracy<T>(
        current_error: T,
        current_grid: T,
        target_error: T,
        order: T,
    ) -> T
    where
        T: RealField + Copy,
    {
        current_grid * (target_error / current_error).powf(T::one() / order)
    }

    /// Compute effective convergence rate between two solutions
    pub fn effective_rate<T>(h1: T, h2: T, error1: T, error2: T) -> T
    where
        T: RealField + Copy,
    {
        (error1 / error2).ln() / (h1 / h2).ln()
    }

    /// Check monotonic convergence
    pub fn is_monotonic<T>(errors: &[T]) -> bool
    where
        T: RealField + Copy,
    {
        errors.windows(2).all(|w| w[1] < w[0])
    }

    /// Compute convergence ratio for oscillatory convergence detection
    pub fn convergence_ratio<T>(errors: &[T]) -> Vec<T>
    where
        T: RealField + Copy,
    {
        if errors.len() < 3 {
            return vec![];
        }

        errors
            .windows(3)
            .map(|w| (w[2] - w[1]) / (w[1] - w[0]))
            .collect()
    }

    /// Detect oscillatory convergence
    pub fn is_oscillatory<T>(errors: &[T]) -> bool
    where
        T: RealField + Copy,
    {
        let ratios = Self::convergence_ratio(errors);

        // Check for sign changes in convergence ratio
        ratios.windows(2).any(|w| w[0] * w[1] < T::zero())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_convergence_order_classification() {
        let order1 = ConvergenceOrder::<f64>::from_rate(0.98);
        assert_eq!(order1, ConvergenceOrder::FirstOrder);

        let order2 = ConvergenceOrder::<f64>::from_rate(2.05);
        assert_eq!(order2, ConvergenceOrder::SecondOrder);

        let order_custom = ConvergenceOrder::<f64>::from_rate(2.5);
        assert!(matches!(order_custom, ConvergenceOrder::Custom(_)));
    }

    #[test]
    fn test_uniform_refinement() {
        let uniform_grids = vec![0.4, 0.2, 0.1, 0.05];
        assert!(ConvergenceAnalysis::is_uniform_refinement(&uniform_grids));

        let non_uniform = vec![0.4, 0.25, 0.1, 0.05];
        assert!(!ConvergenceAnalysis::is_uniform_refinement(&non_uniform));
    }

    #[test]
    fn test_monotonic_convergence() {
        let monotonic = vec![0.1, 0.05, 0.025, 0.01];
        assert!(ConvergenceAnalysis::is_monotonic(&monotonic));

        let oscillatory = vec![0.1, 0.05, 0.06, 0.03];
        assert!(!ConvergenceAnalysis::is_monotonic(&oscillatory));
    }
}
