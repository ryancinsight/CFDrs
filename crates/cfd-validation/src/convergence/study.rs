//! Grid convergence study implementation
//!
//! Implements convergence analysis following Richardson (1911) and Roache (1998) methodologies.

use crate::error_metrics::ErrorAnalysis;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Grid convergence study results
#[derive(Debug, Clone)]
pub struct ConvergenceStudy<T: RealField + Copy> {
    /// Grid sizes used in the study
    pub grid_sizes: Vec<T>,
    /// Corresponding errors or solution values
    pub errors: Vec<T>,
    /// Computed convergence rate (order of accuracy)
    pub convergence_rate: T,
    /// Error coefficient in Richardson extrapolation
    pub error_coefficient: T,
    /// R-squared value for logarithmic fit quality
    pub r_squared: T,
    /// Extrapolated value at zero grid spacing
    pub extrapolated_value: Option<T>,
}

impl<T: RealField + Copy + FromPrimitive> ConvergenceStudy<T> {
    /// Create a new convergence study from grid sizes and errors
    ///
    /// # Arguments
    /// * `grid_sizes` - Characteristic grid sizes (h)
    /// * `errors` - Corresponding errors or solution values
    ///
    /// # Returns
    /// Convergence study with computed convergence rate
    pub fn new(grid_sizes: Vec<T>, errors: Vec<T>) -> Result<Self> {
        if grid_sizes.len() != errors.len() || grid_sizes.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 grid sizes with corresponding errors".to_string()
            ));
        }

        // Verify monotonic grid refinement
        for i in 1..grid_sizes.len() {
            if grid_sizes[i] >= grid_sizes[i - 1] {
                return Err(Error::InvalidConfiguration(
                    "Grid sizes must be monotonically decreasing".to_string()
                ));
            }
        }

        // Compute convergence rate using least squares fit
        let convergence_rate = compute_convergence_rate(&grid_sizes, &errors)?;

        // Compute error coefficient and fit quality
        let (error_coefficient, r_squared) = compute_fit_quality(
            &grid_sizes,
            &errors,
            convergence_rate
        )?;

        Ok(Self {
            grid_sizes,
            errors,
            convergence_rate,
            error_coefficient,
            r_squared,
            extrapolated_value: None,
        })
    }

    /// Check if the study is in the asymptotic range
    ///
    /// Returns true if R² > 0.99, indicating consistent convergence behavior
    pub fn is_asymptotic(&self) -> bool {
        self.r_squared > T::from_f64(0.99).unwrap_or_else(T::one)
    }

    /// Predict error for a given grid size using the power law model
    ///
    /// error = C * h^p where C is error_coefficient and p is convergence_rate
    pub fn predict_error(&self, grid_size: T) -> T {
        self.error_coefficient * grid_size.powf(self.convergence_rate)
    }

    /// Estimate grid size needed to achieve target error
    pub fn grid_size_for_error(&self, target_error: T) -> Result<T> {
        if self.error_coefficient <= T::zero() {
            return Err(Error::InvalidInput(
                "Error coefficient must be positive".to_string()
            ));
        }
        
        if target_error <= T::zero() {
            return Err(Error::InvalidInput(
                "Target error must be positive".to_string()
            ));
        }

        Ok((target_error / self.error_coefficient).powf(T::one() / self.convergence_rate))
    }
}

/// Compute convergence rate using least squares regression
///
/// Fits log(error) = log(C) + p * log(h) to determine convergence order p
pub fn compute_convergence_rate<T>(grid_sizes: &[T], errors: &[T]) -> Result<T>
where
    T: RealField + Copy + FromPrimitive,
{
    let n = T::from_usize(grid_sizes.len())
        .ok_or_else(|| Error::InvalidInput("Cannot convert size".to_string()))?;

    // Compute logarithms for linear regression
    let (sum_log_h, sum_log_e, sum_log_h2, sum_log_he) = 
        grid_sizes.iter()
            .zip(errors.iter())
            .map(|(h, e)| {
                let log_h = h.ln();
                let log_e = e.ln();
                (log_h, log_e, log_h * log_h, log_h * log_e)
            })
            .fold(
                (T::zero(), T::zero(), T::zero(), T::zero()),
                |(sh, se, sh2, she), (lh, le, lh2, lhe)| {
                    (sh + lh, se + le, sh2 + lh2, she + lhe)
                }
            );

    // Least squares solution for slope (convergence rate)
    let denominator = n * sum_log_h2 - sum_log_h * sum_log_h;
    
    if denominator.abs() < T::from_f64(cfd_core::constants::numerical::solver::EPSILON_TOLERANCE)
        .unwrap_or_else(|| T::from_f64(1e-10).unwrap()) 
    {
        return Err(Error::Numerical(
            cfd_core::error::NumericalErrorKind::SingularMatrix
        ));
    }

    let convergence_rate = (n * sum_log_he - sum_log_h * sum_log_e) / denominator;
    
    Ok(convergence_rate)
}

/// Compute fit quality metrics (error coefficient and R-squared)
fn compute_fit_quality<T>(
    grid_sizes: &[T],
    errors: &[T],
    convergence_rate: T,
) -> Result<(T, T)>
where
    T: RealField + Copy + FromPrimitive,
{
    let n = T::from_usize(grid_sizes.len())
        .ok_or_else(|| Error::InvalidInput("Cannot convert size".to_string()))?;

    // Compute error coefficient from intercept
    let (sum_log_h, sum_log_e) = grid_sizes.iter()
        .zip(errors.iter())
        .map(|(h, e)| (h.ln(), e.ln()))
        .fold(
            (T::zero(), T::zero()),
            |(sh, se), (lh, le)| (sh + lh, se + le)
        );

    let log_c = (sum_log_e - convergence_rate * sum_log_h) / n;
    let error_coefficient = log_c.exp();

    // Compute R-squared for fit quality
    let mean_log_e = sum_log_e / n;
    
    let (ss_tot, ss_res) = errors.iter()
        .zip(grid_sizes.iter())
        .map(|(e, h)| {
            let log_e = e.ln();
            let log_h = h.ln();
            let predicted = log_c + convergence_rate * log_h;
            let total_dev = log_e - mean_log_e;
            let residual = log_e - predicted;
            (total_dev * total_dev, residual * residual)
        })
        .fold(
            (T::zero(), T::zero()),
            |(tot, res), (t, r)| (tot + t, res + r)
        );

    let r_squared = if ss_tot > T::zero() {
        T::one() - ss_res / ss_tot
    } else {
        T::one()
    };

    Ok((error_coefficient, r_squared))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_convergence_study() {
        // Test with known second-order convergence
        let grid_sizes = vec![0.1, 0.05, 0.025];
        let errors = vec![0.01, 0.0025, 0.000625]; // O(h²) behavior
        
        let study = ConvergenceStudy::<f64>::new(grid_sizes, errors).unwrap();
        
        assert_relative_eq!(study.convergence_rate, 2.0, epsilon = 0.01);
        assert!(study.r_squared > 0.99);
        assert!(study.is_asymptotic());
    }

    #[test]
    fn test_grid_refinement_ratio() {
        let grid_sizes = vec![0.2, 0.1, 0.05];
        let errors = vec![0.04, 0.01, 0.0025];
        
        let study = ConvergenceStudy::<f64>::new(grid_sizes, errors).unwrap();
        
        // Predict error for h=0.025
        let predicted = study.predict_error(0.025);
        assert_relative_eq!(predicted, 0.000625, epsilon = 1e-6);
    }
}