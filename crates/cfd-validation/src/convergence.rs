//! Convergence analysis tools for CFD validation.
//!
//! This module provides tools for analyzing the convergence behavior of CFD solvers,
//! including Richardson extrapolation and grid convergence studies.

use crate::error_metrics::{ErrorMetric, ErrorAnalysis};
use cfd_core::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::{FromPrimitive, ToPrimitive};

/// Convergence study results
#[derive(Debug, Clone)]
pub struct ConvergenceStudy<T: RealField> {
    /// Grid sizes used in the study
    pub grid_sizes: Vec<T>,
    /// Corresponding errors
    pub errors: Vec<T>,
    /// Computed convergence rate
    pub convergence_rate: T,
    /// Coefficient in the error model: error = C * h^p
    pub error_coefficient: T,
    /// R-squared value for the fit quality
    pub r_squared: T,
    /// Richardson extrapolated value (h → 0)
    pub extrapolated_value: Option<T>,
}

impl<T: RealField + FromPrimitive> ConvergenceStudy<T> {
    /// Create a new convergence study
    pub fn new(grid_sizes: Vec<T>, errors: Vec<T>) -> Result<Self> {
        if grid_sizes.len() != errors.len() || grid_sizes.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 grid sizes and corresponding errors".to_string()
            ));
        }

        // Compute convergence rate
        let convergence_rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors)?;

        // Compute error coefficient and R-squared
        let (error_coefficient, r_squared) = Self::compute_fit_quality(&grid_sizes, &errors, convergence_rate.clone())?;

        Ok(Self {
            grid_sizes,
            errors,
            convergence_rate,
            error_coefficient,
            r_squared,
            extrapolated_value: None,
        })
    }

    /// Compute fit quality metrics using single-pass iterator calculations
    fn compute_fit_quality(
        grid_sizes: &[T],
        errors: &[T],
        convergence_rate: T,
    ) -> Result<(T, T)> {
        let n = T::from_usize(grid_sizes.len()).unwrap();

        // Single-pass calculation of all required statistics
        let (sum_log_e, sum_log_h, sum_log_e_sq, _sum_log_h_sq) = errors.iter()
            .zip(grid_sizes.iter())
            .map(|(e, h)| (e.clone().ln(), h.clone().ln()))
            .fold(
                (T::zero(), T::zero(), T::zero(), T::zero()),
                |(acc_e, acc_h, acc_e_sq, acc_h_sq), (log_e, log_h)| {
                    (
                        acc_e + log_e.clone(),
                        acc_h + log_h.clone(),
                        acc_e_sq + log_e.clone() * log_e,
                        acc_h_sq + log_h.clone() * log_h
                    )
                }
            );

        let log_c = (sum_log_e.clone() - convergence_rate.clone() * sum_log_h) / n.clone();
        let error_coefficient = log_c.clone().exp();

        // R-squared calculation using single-pass statistics
        let mean_log_error = sum_log_e / n.clone();
        let ss_tot = sum_log_e_sq - n.clone() * mean_log_error.clone() * mean_log_error;
        
        // Calculate residual sum of squares in single pass
        let ss_res = errors.iter()
            .zip(grid_sizes.iter())
            .map(|(e, h)| {
                let log_e = e.clone().ln();
                let log_h = h.clone().ln();
                let predicted = log_c.clone() + convergence_rate.clone() * log_h;
                let diff = log_e - predicted;
                diff.clone() * diff
            })
            .fold(T::zero(), |acc, x| acc + x);

        let r_squared = if ss_tot != T::zero() {
            T::one() - ss_res / ss_tot
        } else {
            T::one()
        };

        Ok((error_coefficient, r_squared))
    }

    /// Check if convergence is achieved
    pub fn is_converged(&self, tolerance: T) -> bool {
        if self.errors.len() < 2 {
            return false;
        }

        let latest_error = &self.errors[self.errors.len() - 1];
        *latest_error <= tolerance
    }

    /// Get the convergence order (theoretical vs observed) with default tolerance
    pub fn convergence_order(&self) -> ConvergenceOrder<T> {
        let default_tolerance = T::from_f64(0.1).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        self.convergence_order_with_tolerance(default_tolerance)
    }

    /// Get the convergence order with configurable tolerance
    pub fn convergence_order_with_tolerance(&self, tolerance: T) -> ConvergenceOrder<T> {
        let rate = self.convergence_rate.clone();
        let one = T::one();
        let two = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        let three = T::from_f64(3.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        let four = T::from_f64(4.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;

        if (rate.clone() - one).abs() < tolerance {
            ConvergenceOrder::FirstOrder
        } else if (rate.clone() - two).abs() < tolerance {
            ConvergenceOrder::SecondOrder
        } else if (rate.clone() - three).abs() < tolerance {
            ConvergenceOrder::ThirdOrder
        } else if (rate.clone() - four).abs() < tolerance {
            ConvergenceOrder::FourthOrder
        } else {
            ConvergenceOrder::Other(rate)
        }
    }

    /// Predict error for a given grid size
    pub fn predict_error(&self, grid_size: T) -> T {
        self.error_coefficient.clone() * grid_size.powf(self.convergence_rate.clone())
    }

    /// Estimate grid size needed to achieve target error
    pub fn grid_size_for_error(&self, target_error: T) -> Result<T> {
        if self.error_coefficient == T::zero() {
            return Err(Error::InvalidInput("Error coefficient is zero".to_string()));
        }

        let epsilon = T::from_f64(1e-9).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        if self.convergence_rate.clone().abs() < epsilon {
            return Err(Error::InvalidInput("Convergence rate is zero or near-zero".to_string()));
        }

        if target_error <= T::zero() {
            return Err(Error::InvalidInput("Target error must be positive".to_string()));
        }

        Ok((target_error / self.error_coefficient.clone()).powf(T::one() / self.convergence_rate.clone()))
    }
}

/// Convergence order classification
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceOrder<T: RealField> {
    /// First-order convergence (p ≈ 1)
    FirstOrder,
    /// Second-order convergence (p ≈ 2)
    SecondOrder,
    /// Third-order convergence (p ≈ 3)
    ThirdOrder,
    /// Fourth-order convergence (p ≈ 4)
    FourthOrder,
    /// Other convergence rate
    Other(T),
}

/// Richardson extrapolation for improved accuracy
pub struct RichardsonExtrapolation<T: RealField> {
    /// Convergence rate (order of accuracy)
    pub convergence_rate: T,
}

impl<T: RealField + FromPrimitive + ToPrimitive> RichardsonExtrapolation<T> {
    /// Create new Richardson extrapolation
    pub fn new(convergence_rate: T) -> Self {
        Self { convergence_rate }
    }

    /// Create Richardson extrapolation with second-order accuracy
    pub fn second_order() -> Self {
        Self::new(T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?)
    }

    /// Extrapolate using two solutions on different grids
    pub fn extrapolate_two_grids(
        &self,
        coarse_solution: T,
        fine_solution: T,
        grid_ratio: T, // h_coarse / h_fine
    ) -> T {
        let r_p = grid_ratio.powf(self.convergence_rate.clone());
        (r_p.clone() * fine_solution - coarse_solution) / (r_p - T::one())
    }

    /// Extrapolate using two solutions (backward compatibility)
    pub fn extrapolate(&self, solution1: T, solution2: T) -> T {
        // Assume default grid ratio of 2.0
        let grid_ratio = T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        self.extrapolate_two_grids(solution1, solution2, grid_ratio)
    }

    /// Extrapolate using three solutions (more robust)
    /// Extrapolate using three grids with default tolerance (10%)
    pub fn extrapolate_three_grids(
        &self,
        coarse_solution: T,
        medium_solution: T,
        fine_solution: T,
        grid_ratio: T, // constant ratio between grids
    ) -> Result<T> {
        let default_tolerance = T::from_f64(0.1).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        self.extrapolate_three_grids_with_tolerance(
            coarse_solution, medium_solution, fine_solution, grid_ratio, default_tolerance
        )
    }

    /// Extrapolate using three grids with configurable tolerance
    pub fn extrapolate_three_grids_with_tolerance(
        &self,
        coarse_solution: T,
        medium_solution: T,
        fine_solution: T,
        grid_ratio: T, // constant ratio between grids
        tolerance: T,
    ) -> Result<T> {
        // Use the finest two grids for extrapolation
        let extrapolated = self.extrapolate_two_grids(medium_solution.clone(), fine_solution, grid_ratio.clone());

        // Verify consistency with the coarse grid
        let expected_coarse = self.extrapolate_two_grids(extrapolated.clone(), medium_solution, grid_ratio);
        let relative_error = (expected_coarse - coarse_solution.clone()).abs() / coarse_solution.abs();

        if relative_error > tolerance {
            return Err(Error::ConvergenceFailure(
                "Richardson extrapolation shows poor consistency across grids".to_string()
            ));
        }

        Ok(extrapolated)
    }

    /// Estimate the exact solution using Richardson extrapolation
    pub fn estimate_exact_solution(&self, solutions: &[T], grid_sizes: &[T]) -> Result<T> {
        if solutions.len() != grid_sizes.len() || solutions.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 solutions and corresponding grid sizes".to_string()
            ));
        }

        if solutions.len() == 2 {
            let grid_ratio = grid_sizes[0].clone() / grid_sizes[1].clone();
            Ok(self.extrapolate_two_grids(solutions[0].clone(), solutions[1].clone(), grid_ratio))
        } else {
            // Use the finest three grids
            let n = solutions.len();
            let grid_ratio = grid_sizes[n-2].clone() / grid_sizes[n-1].clone();
            self.extrapolate_three_grids(
                solutions[n-3].clone(),
                solutions[n-2].clone(),
                solutions[n-1].clone(),
                grid_ratio,
            )
        }
    }
}

/// Grid convergence index (GCI) analysis
pub struct GridConvergenceIndex<T: RealField> {
    /// Safety factor (typically 1.25 for 3+ grids, 3.0 for 2 grids)
    pub safety_factor: T,
}

impl<T: RealField + FromPrimitive> GridConvergenceIndex<T> {
    /// Create new GCI analysis
    pub fn new(safety_factor: T) -> Self {
        Self { safety_factor }
    }

    /// Create GCI for three or more grids
    pub fn for_multiple_grids() -> Self {
        Self::new(T::from_f64(1.25).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?)
    }

    /// Create GCI for two grids
    pub fn for_two_grids() -> Self {
        Self::new(T::from_f64(3.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?)
    }

    /// Compute GCI between two solutions
    pub fn compute_gci(
        &self,
        coarse_solution: T,
        fine_solution: T,
        grid_ratio: T,
        convergence_rate: T,
    ) -> T {
        let relative_error = (coarse_solution.clone() - fine_solution.clone()).abs() / fine_solution.abs();
        let r_p = grid_ratio.powf(convergence_rate);

        self.safety_factor.clone() * relative_error / (r_p - T::one())
    }

    /// Check if solutions are in asymptotic range
    pub fn is_asymptotic_range(
        &self,
        gci_coarse: T,
        gci_fine: T,
        grid_ratio: T,
        convergence_rate: T,
    ) -> bool {
        let r_p = grid_ratio.powf(convergence_rate);
        let ratio = gci_coarse / (r_p * gci_fine);

        // Should be close to 1.0 in asymptotic range
        let tolerance = T::from_f64(0.01).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?; // 1% tolerance
        (ratio - T::one()).abs() < tolerance
    }
}

/// Convergence analysis utilities
pub struct ConvergenceAnalysis;

impl ConvergenceAnalysis {
    /// Perform a complete grid convergence study
    pub fn grid_convergence_study<T, M>(
        solutions: &[T],
        grid_sizes: &[T],
        metric: &M,
        reference_solution: Option<&[T]>,
    ) -> Result<ConvergenceStudy<T>>
    where
        T: RealField + FromPrimitive + ToPrimitive + Clone,
        M: ErrorMetric<T>,
    {
        if solutions.len() != grid_sizes.len() || solutions.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 solutions and corresponding grid sizes".to_string()
            ));
        }

        let errors = if let Some(reference) = reference_solution {
            // Compute errors against reference solution
            solutions.iter()
                .map(|sol| metric.compute_error(&[sol.clone()], reference).unwrap_or(T::zero()))
                .collect()
        } else {
            // Use Richardson extrapolation to estimate errors
            let richardson = RichardsonExtrapolation::second_order();
            let mut errors = Vec::new();

            for i in 0..solutions.len() - 1 {
                let grid_ratio = grid_sizes[i].clone() / grid_sizes[i + 1].clone();
                let extrapolated = richardson.extrapolate_two_grids(
                    solutions[i].clone(),
                    solutions[i + 1].clone(),
                    grid_ratio,
                );
                let error = (solutions[i].clone() - extrapolated).abs();
                errors.push(error);
            }

                         // Estimate error for finest grid using Richardson extrapolation
             // For the finest grid (h1) and the next one (h2), the error can be estimated as:
             // error1 ≈ |f1 - f2| / ( (h2/h1)^p - 1 )
             if errors.len() >= 2 && solutions.len() >= 2 {
                 let h1 = grid_sizes[solutions.len() - 1].clone();
                 let h2 = grid_sizes[solutions.len() - 2].clone();
                 let f1 = solutions[solutions.len() - 1].clone();
                 let f2 = solutions[solutions.len() - 2].clone();
                 
                 // Use computed convergence rate from existing data
                 let p = if errors.len() >= 2 {
                     // Estimate convergence rate from last two error points
                     let h_ratio = grid_sizes[errors.len() - 2].clone() / grid_sizes[errors.len() - 1].clone();
                     let e_ratio = errors[errors.len() - 2].clone() / errors[errors.len() - 1].clone();
                     e_ratio.ln() / h_ratio.ln()
                 } else {
                     T::from_f64(2.0).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? // Assume second-order as default
                 };
                 
                 let h_ratio = h2 / h1;
                 let denominator = h_ratio.powf(p) - T::one();
                 
                 if denominator.clone().abs() > T::from_f64(1e-10).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? {
                     let finest_error = (f1 - f2).abs() / denominator;
                     errors.push(finest_error);
                 } else {
                     // Fallback to simple extrapolation if denominator is too small
                     let ratio = errors[errors.len() - 1].clone() / errors[errors.len() - 2].clone();
                     errors.push(errors[errors.len() - 1].clone() * ratio);
                 }
             } else {
                 // Simple fallback for insufficient data
                 errors.push(errors[0].clone() * T::from_f64(0.5).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?);
             }

            errors
        };

        ConvergenceStudy::new(grid_sizes.to_vec(), errors)
    }

    /// Recommend next grid size for convergence study
    pub fn recommend_next_grid_size<T: RealField + FromPrimitive>(
        current_grid_sizes: &[T],
        target_error_reduction: T,
        convergence_rate: T,
    ) -> T {
        if current_grid_sizes.is_empty() {
            return T::one();
        }

        let finest_grid = current_grid_sizes[current_grid_sizes.len() - 1].clone();
        let refinement_ratio = target_error_reduction.powf(T::one() / convergence_rate);

        finest_grid / refinement_ratio
    }

    /// Check convergence criteria
    pub fn check_convergence<T: RealField>(
        errors: &[T],
        absolute_tolerance: T,
        relative_tolerance: T,
    ) -> ConvergenceStatus<T> {
        if errors.is_empty() {
            return ConvergenceStatus::InsufficientData;
        }

        let latest_error = &errors[errors.len() - 1];

        // Check absolute convergence
        if *latest_error <= absolute_tolerance {
            return ConvergenceStatus::Converged {
                final_error: latest_error.clone(),
                criterion: ConvergenceCriterion::Absolute,
            };
        }

        // Check relative convergence (if we have multiple errors)
        if errors.len() >= 2 {
            let previous_error = &errors[errors.len() - 2];
            let relative_change = (latest_error.clone() - previous_error.clone()).abs() / previous_error.abs();

            if relative_change <= relative_tolerance {
                return ConvergenceStatus::Converged {
                    final_error: latest_error.clone(),
                    criterion: ConvergenceCriterion::Relative,
                };
            }
        }

        ConvergenceStatus::NotConverged {
            current_error: latest_error.clone(),
        }
    }
}

/// Convergence status
#[derive(Debug, Clone, PartialEq)]
pub enum ConvergenceStatus<T: RealField> {
    /// Convergence achieved
    Converged {
        /// Final error value
        final_error: T,
        /// Convergence criterion that was satisfied
        criterion: ConvergenceCriterion,
    },
    /// Convergence not yet achieved
    NotConverged {
        /// Current error value
        current_error: T,
    },
    /// Insufficient data to determine convergence
    InsufficientData,
}

/// Convergence criteria
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ConvergenceCriterion {
    /// Absolute error tolerance
    Absolute,
    /// Relative error tolerance
    Relative,
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_convergence_study_second_order() {
        // Test with theoretical second-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625]; // errors ∝ h^2

        let study = ConvergenceStudy::new(grid_sizes, errors).unwrap();

        assert_relative_eq!(study.convergence_rate, 2.0, epsilon = 1e-10);
        assert!(study.r_squared > 0.99); // Should have excellent fit

        match study.convergence_order() {
            ConvergenceOrder::SecondOrder => {},
            _ => assert!(false, "Expected second-order convergence"),
        }
    }

    #[test]
    fn test_convergence_study_first_order() {
        // Test with theoretical first-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.1, 0.05, 0.025]; // errors ∝ h^1

        let study = ConvergenceStudy::new(grid_sizes, errors).unwrap();

        assert_relative_eq!(study.convergence_rate, 1.0, epsilon = 1e-10);

        match study.convergence_order() {
            ConvergenceOrder::FirstOrder => {},
            _ => assert!(false, "Expected first-order convergence"),
        }
    }

    #[test]
    fn test_error_prediction() {
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625];

        let study = ConvergenceStudy::new(grid_sizes, errors).unwrap();

        // Predict error for h = 0.0125 (should be ~0.000156)
        let predicted = study.predict_error(0.0125);
        let expected = 0.01_f64 * (0.0125_f64 / 0.1).powi(2); // C * h^2
        assert_relative_eq!(predicted, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_grid_size_for_error() {
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625];

        let study = ConvergenceStudy::new(grid_sizes, errors).unwrap();

        // Find grid size for target error of 0.0001
        let target_error = 0.0001;
        let required_h = study.grid_size_for_error(target_error).unwrap();

        // Verify by predicting error at this grid size
        let predicted_error = study.predict_error(required_h);
        assert_relative_eq!(predicted_error, target_error, epsilon = 1e-6);
    }

    #[test]
    fn test_grid_size_for_error_edge_cases() {
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625];
        let study = ConvergenceStudy::new(grid_sizes, errors).unwrap();

        // Test with zero target error
        assert!(study.grid_size_for_error(0.0).is_err());
        
        // Test with negative target error
        assert!(study.grid_size_for_error(-0.001).is_err());
    }

    #[test]
    fn test_richardson_extrapolation_two_grids() {
        let richardson = RichardsonExtrapolation::second_order();

        // Test with known second-order convergence
        let coarse_solution = 1.01; // exact + 0.01 error
        let fine_solution = 1.0025; // exact + 0.0025 error
        let grid_ratio = 2.0; // h_coarse / h_fine

        let extrapolated = richardson.extrapolate_two_grids(coarse_solution, fine_solution, grid_ratio);

        // Should be close to exact solution (1.0)
        assert_relative_eq!(extrapolated, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_richardson_extrapolation_three_grids() {
        let richardson = RichardsonExtrapolation::second_order();

        let coarse_solution = 1.04;   // exact + 0.04 error
        let medium_solution = 1.01;   // exact + 0.01 error
        let fine_solution = 1.0025;   // exact + 0.0025 error
        let grid_ratio = 2.0;

        let extrapolated = richardson.extrapolate_three_grids(
            coarse_solution,
            medium_solution,
            fine_solution,
            grid_ratio,
        ).unwrap();

        // Should be close to exact solution (1.0)
        assert_relative_eq!(extrapolated, 1.0, epsilon = 1e-6);
    }

    #[test]
    fn test_gci_computation() {
        let gci = GridConvergenceIndex::for_multiple_grids();

        let coarse_solution = 1.01;
        let fine_solution = 1.0025;
        let grid_ratio = 2.0;
        let convergence_rate = 2.0;

        let gci_value = gci.compute_gci(coarse_solution, fine_solution, grid_ratio, convergence_rate);

        assert!(gci_value > 0.0);
        assert!(gci_value < 1.0); // Should be reasonable
    }

    #[test]
    fn test_convergence_status_absolute() {
        let errors = vec![0.1, 0.05, 0.001];
        let status = ConvergenceAnalysis::check_convergence(&errors, 0.01, 0.1);

        match status {
            ConvergenceStatus::Converged { criterion: ConvergenceCriterion::Absolute, .. } => {},
            _ => assert!(false, "Expected absolute convergence"),
        }
    }

    #[test]
    fn test_convergence_status_relative() {
        let errors = vec![0.1, 0.05, 0.049];
        let status = ConvergenceAnalysis::check_convergence(&errors, 0.001, 0.1);

        match status {
            ConvergenceStatus::Converged { criterion: ConvergenceCriterion::Relative, .. } => {},
            _ => assert!(false, "Expected relative convergence"),
        }
    }

    #[test]
    fn test_convergence_status_not_converged() {
        let errors = vec![0.1, 0.05, 0.04];
        let status = ConvergenceAnalysis::check_convergence(&errors, 0.001, 0.01);

        match status {
            ConvergenceStatus::NotConverged { .. } => {},
            _ => assert!(false, "Expected not converged"),
        }
    }

    #[test]
    fn test_recommend_next_grid_size() {
        let current_grids = vec![0.1, 0.05, 0.025];
        let target_reduction = 4.0; // Want 4x error reduction
        let convergence_rate = 2.0;

        let next_grid = ConvergenceAnalysis::recommend_next_grid_size(
            &current_grids,
            target_reduction,
            convergence_rate,
        );

        // For second-order convergence, 4x error reduction needs 2x grid refinement
        let expected = 0.025 / 2.0;
        assert_relative_eq!(next_grid, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_insufficient_data() {
        let grid_sizes = vec![0.1];
        let errors = vec![0.01];

        assert!(ConvergenceStudy::new(grid_sizes, errors).is_err());
    }

    #[test]
    fn test_mismatched_lengths() {
        let grid_sizes = vec![0.1, 0.05];
        let errors = vec![0.01, 0.0025, 0.000625];

        assert!(ConvergenceStudy::new(grid_sizes, errors).is_err());
    }
}
