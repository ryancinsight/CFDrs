//! Quality metric definitions.

use crate::core::scalar::Real;

/// A quality metric measurement.
#[derive(Clone, Debug)]
pub struct QualityMetric {
    /// Minimum value across all elements.
    pub min: Real,
    /// Maximum value across all elements.
    pub max: Real,
    /// Mean value across all elements.
    pub mean: Real,
    /// Number of elements measured.
    pub count: usize,
}

impl QualityMetric {
    /// Create from a slice of values.
    pub fn from_values(values: &[Real]) -> Option<Self> {
        if values.is_empty() {
            return None;
        }
        let min = values.iter().copied().fold(Real::INFINITY, Real::min);
        let max = values.iter().copied().fold(Real::NEG_INFINITY, Real::max);
        let sum: Real = values.iter().sum();
        Some(Self {
            min,
            max,
            mean: sum / values.len() as Real,
            count: values.len(),
        })
    }

    /// Number of elements below a threshold.
    pub fn count_below(values: &[Real], threshold: Real) -> usize {
        values.iter().filter(|&&v| v < threshold).count()
    }
}
