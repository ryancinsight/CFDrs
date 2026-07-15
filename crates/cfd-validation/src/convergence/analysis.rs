//! Convergence analysis utilities
//!
//! Provides tools for analyzing and classifying convergence behavior.

use crate::scalar;
use eunomia::{FloatElement, RealField};

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

impl<T: RealField + Copy + FloatElement> ConvergenceOrder<T> {
    /// Determine convergence order from observed rate
    pub fn from_rate(rate: T) -> Self {
        let tolerance = scalar::from_f64::<T>(0.15);

        if rate < scalar::from_f64::<T>(0.5) {
            Self::SubLinear
        } else if scalar::abs(rate - scalar::one::<T>()) < tolerance {
            Self::FirstOrder
        } else if scalar::abs(rate - scalar::from_f64::<T>(2.0)) < tolerance {
            Self::SecondOrder
        } else if scalar::abs(rate - scalar::from_f64::<T>(3.0)) < tolerance {
            Self::ThirdOrder
        } else if scalar::abs(rate - scalar::from_f64::<T>(4.0)) < tolerance {
            Self::FourthOrder
        } else if rate > scalar::from_f64::<T>(6.0) {
            Self::Spectral
        } else if rate > scalar::one::<T>() && rate < scalar::from_f64::<T>(2.0) {
            Self::SuperLinear
        } else {
            Self::Custom(rate)
        }
    }

    /// Get expected rate for this order
    pub fn expected_rate(&self) -> T {
        match self {
            Self::SubLinear => scalar::from_f64::<T>(0.5),
            Self::FirstOrder => scalar::one::<T>(),
            Self::SuperLinear => scalar::from_f64::<T>(1.5),
            Self::SecondOrder => scalar::from_f64::<T>(2.0),
            Self::ThirdOrder => scalar::from_f64::<T>(3.0),
            Self::FourthOrder => scalar::from_f64::<T>(4.0),
            Self::Spectral => scalar::from_f64::<T>(10.0),
            Self::Custom(rate) => *rate,
        }
    }

    /// Check if observed rate matches expected within tolerance
    pub fn matches(&self, observed_rate: T) -> bool {
        let expected = self.expected_rate();
        let tolerance = scalar::from_f64::<T>(0.25); // 25% tolerance
        scalar::abs(observed_rate - expected) < expected * tolerance
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
        T: RealField + Copy + FloatElement + std::iter::Sum,
    {
        let ratios = Self::refinement_ratio(grid_sizes);
        if ratios.is_empty() {
            return true;
        }

        let mean_ratio = ratios.iter().copied().sum::<T>() / scalar::from_usize::<T>(ratios.len());
        let tolerance = scalar::from_f64::<T>(0.05); // 5% variation allowed

        ratios
            .iter()
            .all(|r| scalar::abs((*r - mean_ratio) / mean_ratio) < tolerance)
    }

    /// Estimate required grid size for target accuracy
    pub fn estimate_grid_for_accuracy<T>(
        current_error: T,
        current_grid: T,
        target_error: T,
        order: T,
    ) -> T
    where
        T: RealField + Copy + FloatElement,
    {
        current_grid * scalar::powf(target_error / current_error, scalar::one::<T>() / order)
    }

    /// Compute effective convergence rate between two solutions
    pub fn effective_rate<T>(h1: T, h2: T, error1: T, error2: T) -> T
    where
        T: RealField + Copy + FloatElement,
    {
        scalar::ln(error1 / error2) / scalar::ln(h1 / h2)
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
        ratios.windows(2).any(|w| w[0] * w[1] < scalar::zero::<T>())
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
