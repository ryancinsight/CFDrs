//! Statistical analysis for mesh quality metrics

use nalgebra::RealField;
use num_traits::Float;
use serde::{Deserialize, Serialize};

/// Statistical summary of quality metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityStatistics<T: RealField> {
    /// Minimum value
    pub min: T,
    /// Maximum value  
    pub max: T,
    /// Mean value
    pub mean: T,
    /// Standard deviation
    pub std_dev: T,
    /// Median value
    pub median: T,
    /// 25th percentile
    pub q1: T,
    /// 75th percentile
    pub q3: T,
    /// Number of samples
    pub count: usize,
}

impl<T: RealField + Float> QualityStatistics<T> {
    /// Create new statistics from samples
    pub fn from_samples(samples: &[T]) -> Self {
        if samples.is_empty() {
            return Self::default();
        }
        
        let mut sorted = samples.to_vec();
        sorted.sort_by(|a, b| a.partial_cmp(b).expect("NaN in samples"));
        
        let count = samples.len();
        let sum: T = samples.iter().cloned().sum();
        let mean = sum / T::from(count).unwrap_or_else(|_| T::one());
        
        let variance: T = samples.iter()
            .map(|x| (*x - mean) * (*x - mean))
            .sum::<T>() / T::from(count - 1).unwrap_or_else(|_| T::one());
        
        let std_dev = variance.sqrt();
        
        let median = if count % 2 == 0 {
            (sorted[count / 2 - 1] + sorted[count / 2]) / (T::one() + T::one())
        } else {
            sorted[count / 2]
        };
        
        let q1_idx = count / 4;
        let q3_idx = 3 * count / 4;
        
        Self {
            min: sorted[0],
            max: sorted[count - 1],
            mean,
            std_dev,
            median,
            q1: sorted[q1_idx],
            q3: sorted[q3_idx],
            count,
        }
    }
}

impl<T: RealField> Default for QualityStatistics<T> {
    fn default() -> Self {
        Self {
            min: T::zero(),
            max: T::zero(),
            mean: T::zero(),
            std_dev: T::zero(),
            median: T::zero(),
            q1: T::zero(),
            q3: T::zero(),
            count: 0,
        }
    }
}

/// Running statistics calculator (Welford's algorithm)
pub struct RunningStats<T: RealField> {
    count: usize,
    mean: T,
    m2: T,
    min: T,
    max: T,
}

impl<T: RealField + Float> RunningStats<T> {
    /// Create new running statistics
    pub fn new() -> Self {
        Self {
            count: 0,
            mean: T::zero(),
            m2: T::zero(),
            min: T::infinity(),
            max: T::neg_infinity(),
        }
    }
    
    /// Add a sample
    pub fn push(&mut self, value: T) {
        self.count += 1;
        let delta = value - self.mean;
        self.mean = self.mean + delta / T::from(self.count).unwrap_or_else(|_| T::one());
        let delta2 = value - self.mean;
        self.m2 = self.m2 + delta * delta2;
        
        self.min = self.min.min(value);
        self.max = self.max.max(value);
    }
    
    /// Get current statistics
    pub fn statistics(&self) -> QualityStatistics<T> {
        let variance = if self.count > 1 {
            self.m2 / T::from(self.count - 1).unwrap_or_else(|_| T::one())
        } else {
            T::zero()
        };
        
        QualityStatistics {
            min: self.min,
            max: self.max,
            mean: self.mean,
            std_dev: variance.sqrt(),
            median: self.mean, // Approximation
            q1: self.mean - variance.sqrt(), // Approximation
            q3: self.mean + variance.sqrt(), // Approximation
            count: self.count,
        }
    }
}