//! Standard norm-based error metrics (L1, L2, L∞)

use super::ErrorMetric;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::cast::FromPrimitive;

/// L2 (Euclidean) norm error metric
pub struct L2Norm;

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for L2Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let sum_squared_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| {
                let diff = *num - *ref_val;
                diff * diff
            })
            .fold(T::zero(), |acc, x| acc + x);

        let n = T::from_usize(numerical.len()).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
        })?;
        Ok((sum_squared_diff / n).sqrt())
    }

    fn name(&self) -> &str {
        "L2 Norm"
    }
}

/// L∞ (maximum) norm error metric
pub struct LInfNorm;

impl<T: RealField + Copy> ErrorMetric<T> for LInfNorm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let max_diff = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| (*num - *ref_val).abs())
            .fold(T::zero(), |acc, x| if x > acc { x } else { acc });

        Ok(max_diff)
    }

    fn name(&self) -> &str {
        "L∞ Norm"
    }
}

/// L1 (Manhattan) norm error metric
pub struct L1Norm;

impl<T: RealField + Copy + FromPrimitive + Copy> ErrorMetric<T> for L1Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string()
            ));
        }

        if numerical.is_empty() {
            return Ok(T::zero());
        }

        let sum_abs_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| (*num - *ref_val).abs())
            .fold(T::zero(), |acc, x| acc + x);

        let n = T::from_usize(numerical.len()).ok_or_else(|| {
            Error::InvalidConfiguration("Failed to convert array length to target type".to_string())
        })?;
        Ok(sum_abs_diff / n)
    }

    fn name(&self) -> &str {
        "L1 Norm"
    }
}