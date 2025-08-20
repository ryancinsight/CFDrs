//! Quality criteria and thresholds for mesh validation

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use super::QualityMetrics;

/// Quality criteria for mesh validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityCriteria<T: RealField + Copy> {
    /// Thresholds for different quality levels
    pub thresholds: QualityThresholds<T>,
    /// Strict mode (fail on any violation)
    pub strict: bool,
}

/// Quality thresholds
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityThresholds<T: RealField> {
    /// Maximum acceptable aspect ratio
    pub max_aspect_ratio: T,
    /// Maximum acceptable skewness
    pub max_skewness: T,
    /// Minimum acceptable orthogonality
    pub min_orthogonality: T,
    /// Minimum acceptable Jacobian
    pub min_jacobian: T,
    /// Minimum overall quality score
    pub min_quality_score: T,
}

impl<T: RealField + Copy> QualityCriteria<T> {
    /// Create criteria with default thresholds
    pub fn default_cfd() -> Self {
        Self {
            thresholds: QualityThresholds::default_cfd(),
            strict: false,
        }
    }
    
    /// Create strict criteria for high-accuracy simulations
    pub fn strict_cfd() -> Self {
        Self {
            thresholds: QualityThresholds::strict_cfd(),
            strict: true,
        }
    }
    
    /// Check if metrics are acceptable
    pub fn is_acceptable(&self, metrics: &QualityMetrics<T>) -> bool {
        let mut acceptable = true;
        
        if metrics.aspect_ratio > self.thresholds.max_aspect_ratio {
            acceptable = false;
        }
        
        if metrics.skewness > self.thresholds.max_skewness {
            acceptable = false;
        }
        
        if metrics.orthogonality < self.thresholds.min_orthogonality {
            acceptable = false;
        }
        
        if metrics.jacobian < self.thresholds.min_jacobian {
            acceptable = false;
        }
        
        if metrics.overall_quality_score < self.thresholds.min_quality_score {
            acceptable = false;
        }
        
        acceptable
    }
}

impl<T: RealField + Copy> QualityThresholds<T> {
    /// Default thresholds for CFD meshes
    pub fn default_cfd() -> Self {
        Self {
            max_aspect_ratio: T::from_f64(10.0).unwrap_or_else(|| T::one()),
            max_skewness: T::from_f64(0.8).unwrap_or_else(|| T::one()),
            min_orthogonality: T::from_f64(0.2).unwrap_or_else(|| T::zero()),
            min_jacobian: T::from_f64(0.1).unwrap_or_else(|| T::zero()),
            min_quality_score: T::from_f64(0.3).unwrap_or_else(|| T::zero()),
        }
    }
    
    /// Strict thresholds for high-accuracy simulations
    pub fn strict_cfd() -> Self {
        Self {
            max_aspect_ratio: T::from_f64(5.0).unwrap_or_else(|| T::one()),
            max_skewness: T::from_f64(0.5).unwrap_or_else(|| T::one()),
            min_orthogonality: T::from_f64(0.5).unwrap_or_else(|| T::zero()),
            min_jacobian: T::from_f64(0.3).unwrap_or_else(|| T::zero()),
            min_quality_score: T::from_f64(0.5).unwrap_or_else(|| T::zero()),
        }
    }
}