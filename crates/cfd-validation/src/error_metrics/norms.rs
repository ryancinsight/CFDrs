//! Standard norm-based error metrics (L1, L2, L∞)

use super::ErrorMetric;
use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, RealField};

/// L2 (Euclidean) norm error metric
pub struct L2Norm;

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for L2Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string(),
            ));
        }

        if numerical.is_empty() {
            return Ok(scalar::zero());
        }

        let sum_squared_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| {
                let diff = *num - *ref_val;
                diff * diff
            })
            .fold(scalar::zero(), |acc, x| acc + x);

        let n = scalar::from_usize(numerical.len());
        Ok(scalar::sqrt(sum_squared_diff / n))
    }

    fn name(&self) -> &'static str {
        "L2 Norm"
    }
}

/// L∞ (maximum) norm error metric
pub struct LInfNorm;

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for LInfNorm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string(),
            ));
        }

        if numerical.is_empty() {
            return Ok(scalar::zero());
        }

        let max_diff = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| scalar::abs(*num - *ref_val))
            .fold(scalar::zero(), |acc, x| scalar::max(acc, x));

        Ok(max_diff)
    }

    fn name(&self) -> &'static str {
        "L∞ Norm"
    }
}

/// L1 (Manhattan) norm error metric
pub struct L1Norm;

impl<T: RealField + FloatElement + Copy> ErrorMetric<T> for L1Norm {
    fn compute_error(&self, numerical: &[T], reference: &[T]) -> Result<T> {
        if numerical.len() != reference.len() {
            return Err(Error::InvalidConfiguration(
                "Arrays must have the same length".to_string(),
            ));
        }

        if numerical.is_empty() {
            return Ok(scalar::zero());
        }

        let sum_abs_diff: T = numerical
            .iter()
            .zip(reference.iter())
            .map(|(num, ref_val)| scalar::abs(*num - *ref_val))
            .fold(scalar::zero(), |acc, x| acc + x);

        let n = scalar::from_usize(numerical.len());
        Ok(sum_abs_diff / n)
    }

    fn name(&self) -> &'static str {
        "L1 Norm"
    }
}
