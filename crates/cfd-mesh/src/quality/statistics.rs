//! Statistical analysis for mesh quality metrics

use cfd_core::error::NumericalErrorKind;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Statistical analysis for mesh quality metrics
pub struct QualityStatistics<T: RealField + Copy> {
    samples: Vec<T>,
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> Default for QualityStatistics<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> QualityStatistics<T> {
    /// Create new statistics collector
    pub fn new() -> Self {
        Self {
            samples: Vec::new(),
        }
    }

    /// Create from samples
    pub fn from_samples(samples: Vec<T>) -> Self {
        Self { samples }
    }

    /// Add a sample to the statistics
    pub fn add_sample(&mut self, value: T) {
        self.samples.push(value);
    }

    /// Get the mean of all samples
    pub fn mean(&self) -> Option<T> {
        if self.samples.is_empty() {
            return None;
        }

        let sum: T = self.samples.iter().copied().sum();
        let count = T::from_usize(self.samples.len())?;
        Some(sum / count)
    }

    /// Get the median of all samples
    pub fn median(&self) -> Result<Option<T>> {
        if self.samples.is_empty() {
            return Ok(None);
        }

        let mut sorted = self.samples.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

        let mid = sorted.len() / 2;
        if sorted.len() % 2 == 0 {
            let a = sorted[mid - 1];
            let b = sorted[mid];
            let two = T::from_f64(2.0).ok_or_else(|| {
                Error::Numerical(NumericalErrorKind::InvalidValue {
                    value: "Cannot convert 2.0".to_string(),
                })
            })?;
            Ok(Some((a + b) / two))
        } else {
            Ok(Some(sorted[mid]))
        }
    }

    /// Get the standard deviation
    pub fn std_dev(&self) -> Result<Option<T>> {
        let mean = match self.mean() {
            Some(m) => m,
            None => return Ok(None),
        };

        if self.samples.len() < 2 {
            return Ok(None);
        }

        let variance: T = self
            .samples
            .iter()
            .map(|&x| {
                let diff = x - mean;
                diff * diff
            })
            .sum();

        let n = T::from_usize(self.samples.len() - 1).ok_or_else(|| {
            Error::Numerical(NumericalErrorKind::InvalidValue {
                value: "Cannot convert sample count".to_string(),
            })
        })?;
        Ok(Some((variance / n).sqrt()))
    }

    /// Get minimum value
    pub fn min(&self) -> Option<T> {
        self.samples
            .iter()
            .copied()
            .min_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Get maximum value
    pub fn max(&self) -> Option<T> {
        self.samples
            .iter()
            .copied()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
    }

    /// Get the number of samples
    pub fn count(&self) -> usize {
        self.samples.len()
    }

    /// Clear all samples
    pub fn clear(&mut self) {
        self.samples.clear();
    }
}
