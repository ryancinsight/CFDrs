//! Core quality metrics for mesh elements
//!
//! Based on established metrics from:
//! - Knupp, P. (2001). "Algebraic mesh quality metrics"
//! - Shewchuk, J. (2002). "What is a good linear finite element?"

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Comprehensive quality metrics for mesh elements
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics<T: RealField> {
    /// Aspect ratio (1.0 = perfect, >1 = stretched)
    pub aspect_ratio: T,
    
    /// Skewness (0.0 = perfect, 1.0 = degenerate)
    pub skewness: T,
    
    /// Orthogonality (1.0 = perfect, 0.0 = non-orthogonal)
    pub orthogonality: T,
    
    /// Smoothness (variation between neighboring cells)
    pub smoothness: T,
    
    /// Jacobian determinant (positive = valid)
    pub jacobian: T,
    
    /// Condition number (1.0 = perfect, higher = worse)
    pub condition_number: T,
    
    /// Volume/Area ratio to ideal element
    pub size_ratio: T,
    
    /// Overall quality score (0-1, higher is better)
    pub overall_quality_score: T,
}

impl<T: RealField> QualityMetrics<T> {
    /// Create metrics with all values set to ideal
    pub fn ideal() -> Self {
        Self {
            aspect_ratio: T::one(),
            skewness: T::zero(),
            orthogonality: T::one(),
            smoothness: T::one(),
            jacobian: T::one(),
            condition_number: T::one(),
            size_ratio: T::one(),
            overall_quality_score: T::one(),
        }
    }
    
    /// Calculate overall quality score from individual metrics
    pub fn calculate_overall_score(&mut self) {
        // Weighted harmonic mean of metrics
        let weights = [
            T::from_f64(0.2).unwrap_or_else(|| T::one()), // aspect_ratio
            T::from_f64(0.2).unwrap_or_else(|| T::one()), // skewness
            T::from_f64(0.15).unwrap_or_else(|| T::one()), // orthogonality
            T::from_f64(0.1).unwrap_or_else(|| T::one()), // smoothness
            T::from_f64(0.15).unwrap_or_else(|| T::one()), // jacobian
            T::from_f64(0.1).unwrap_or_else(|| T::one()), // condition
            T::from_f64(0.1).unwrap_or_else(|| T::one()), // size_ratio
        ];
        
        let metrics = [
            self.aspect_ratio_quality(),
            T::one() - self.skewness,
            self.orthogonality,
            self.smoothness,
            self.jacobian_quality(),
            self.condition_quality(),
            self.size_ratio_quality(),
        ];
        
        let mut weighted_sum = T::zero();
        let mut weight_total = T::zero();
        
        for (metric, weight) in metrics.iter().zip(weights.iter()) {
            if *metric > T::zero() {
                weighted_sum = weighted_sum + *weight / *metric;
                weight_total = weight_total + *weight;
            }
        }
        
        self.overall_quality_score = if weighted_sum > T::zero() {
            weight_total / weighted_sum
        } else {
            T::zero()
        };
    }
    
    /// Convert aspect ratio to quality metric (0-1)
    fn aspect_ratio_quality(&self) -> T {
        T::one() / self.aspect_ratio.max(T::one() / self.aspect_ratio)
    }
    
    /// Convert Jacobian to quality metric (0-1)
    fn jacobian_quality(&self) -> T {
        if self.jacobian > T::zero() {
            (T::one() - (T::one() - self.jacobian).abs()).max(T::zero())
        } else {
            T::zero()
        }
    }
    
    /// Convert condition number to quality metric (0-1)
    fn condition_quality(&self) -> T {
        T::one() / self.condition_number
    }
    
    /// Convert size ratio to quality metric (0-1)
    fn size_ratio_quality(&self) -> T {
        (T::one() - (T::one() - self.size_ratio).abs()).max(T::zero())
    }
}