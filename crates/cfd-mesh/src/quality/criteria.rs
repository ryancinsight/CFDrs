//! Quality criteria and thresholds for mesh validation

use super::QualityMetrics;
use cfd_core::numeric;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
/// Quality criteria for mesh validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityCriteria<T: RealField + Copy> {
    /// Thresholds for different quality levels
    pub thresholds: QualityThresholds<T>,
    /// Strict mode (fail on any violation)
    pub strict: bool,
}
/// Quality thresholds
pub struct QualityThresholds<T: RealField + Copy> {
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
impl<T: RealField + Copy> QualityCriteria<T> {
    /// Create criteria with default thresholds
    #[must_use]
    pub fn default_cfd() -> Self {
        Self {
            thresholds: QualityThresholds::default_cfd(),
            strict: false,
        }
    }
    /// Create strict criteria for high-accuracy simulations
    pub fn strict_cfd() -> Self {
            thresholds: QualityThresholds::strict_cfd(),
            strict: true,
    /// Check if metrics are acceptable
    pub fn is_acceptable(&self, metrics: &QualityMetrics<T>) -> bool {
        let mut acceptable = true;
        if metrics.aspect_ratio > self.thresholds.max_aspect_ratio {
            acceptable = false;
        if metrics.skewness > self.thresholds.max_skewness {
        if metrics.orthogonality < self.thresholds.min_orthogonality {
        if metrics.jacobian < self.thresholds.min_jacobian {
        if metrics.overall_quality_score < self.thresholds.min_quality_score {
        acceptable
impl<T: RealField + Copy> QualityThresholds<T> {
    /// Default thresholds for CFD meshes
            max_aspect_ratio: cfd_core::numeric::from_f64(10.0)?,
            max_skewness: cfd_core::numeric::from_f64(0.8)?,
            min_orthogonality: cfd_core::numeric::from_f64(0.2)?,
            min_jacobian: cfd_core::numeric::from_f64(0.1)?,
            min_quality_score: cfd_core::numeric::from_f64(0.3)?,
    /// Strict thresholds for high-accuracy simulations
            max_aspect_ratio: cfd_core::numeric::from_f64(5.0)?,
            max_skewness: cfd_core::numeric::from_f64(0.5)?,
            min_orthogonality: cfd_core::numeric::from_f64(0.5)?,
            min_jacobian: cfd_core::numeric::from_f64(0.3)?,
            min_quality_score: cfd_core::numeric::from_f64(0.5)?,
