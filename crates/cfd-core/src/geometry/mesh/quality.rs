//! Mesh quality assessment and metrics
//!
//! This module provides tools for evaluating the geometric quality of mesh
//! elements, which is critical for simulation stability and accuracy.

use super::Mesh;
use crate::error::Result;
use eunomia::FloatElement;
use eunomia::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Service for mesh quality assessment
pub struct MeshQualityService;

impl MeshQualityService {
    /// Assess overall mesh quality and provide recommendations
    pub fn assess_quality<T: RealField + FloatElement + Copy>(
        aspect_ratio_stats: &MetricStatistics<T>,
        skewness_stats: &MetricStatistics<T>,
        orthogonality_stats: &MetricStatistics<T>,
    ) -> QualityAssessment {
        let aspect_ratio_quality = Self::assess_aspect_ratio(aspect_ratio_stats);
        let skewness_quality = Self::assess_skewness(skewness_stats);
        let orthogonality_quality = Self::assess_orthogonality(orthogonality_stats);

        let overall_quality = [
            aspect_ratio_quality,
            skewness_quality,
            orthogonality_quality,
        ]
        .iter()
        .min()
        .copied()
        .unwrap_or(QualityLevel::Level1);

        QualityAssessment {
            overall_quality,
            aspect_ratio_quality,
            skewness_quality,
            orthogonality_quality,
            recommendations: Self::generate_recommendations(
                aspect_ratio_quality,
                skewness_quality,
                orthogonality_quality,
            ),
        }
    }

    fn assess_aspect_ratio<T: RealField + FloatElement + Copy>(
        stats: &MetricStatistics<T>,
    ) -> QualityLevel {
        let thresh_level4 = <T as FloatElement>::from_f64(2.0);
        let thresh_level3 = <T as FloatElement>::from_f64(5.0);
        let thresh_level2 = <T as FloatElement>::from_f64(10.0);

        if stats.max < thresh_level4 {
            QualityLevel::Level4
        } else if stats.max < thresh_level3 {
            QualityLevel::Level3
        } else if stats.max < thresh_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_skewness<T: RealField + FloatElement + Copy>(
        stats: &MetricStatistics<T>,
    ) -> QualityLevel {
        let thresh_level4 = <T as FloatElement>::from_f64(0.25);
        let thresh_level3 = <T as FloatElement>::from_f64(0.5);
        let thresh_level2 = <T as FloatElement>::from_f64(0.8);

        if stats.max < thresh_level4 {
            QualityLevel::Level4
        } else if stats.max < thresh_level3 {
            QualityLevel::Level3
        } else if stats.max < thresh_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_orthogonality<T: RealField + FloatElement + Copy>(
        stats: &MetricStatistics<T>,
    ) -> QualityLevel {
        let thresh_level4 = <T as FloatElement>::from_f64(0.95);
        let thresh_level3 = <T as FloatElement>::from_f64(0.85);
        let thresh_level2 = <T as FloatElement>::from_f64(0.7);

        if stats.min > thresh_level4 {
            QualityLevel::Level4
        } else if stats.min > thresh_level3 {
            QualityLevel::Level3
        } else if stats.min > thresh_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn generate_recommendations(
        aspect_ratio: QualityLevel,
        skewness: QualityLevel,
        orthogonality: QualityLevel,
    ) -> Vec<String> {
        let mut recommendations = Vec::new();

        if aspect_ratio == QualityLevel::Level1 {
            recommendations
                .push("Consider refining mesh in regions with high aspect ratio cells".to_string());
        }

        if skewness == QualityLevel::Level1 {
            recommendations.push("Improve mesh quality by reducing cell skewness".to_string());
        }

        if orthogonality == QualityLevel::Level1 {
            recommendations
                .push("Enhance mesh orthogonality for better numerical accuracy".to_string());
        }

        if recommendations.is_empty() {
            recommendations.push("Mesh quality is acceptable for CFD simulation".to_string());
        }

        recommendations
    }
}

/// Quality level enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum QualityLevel {
    /// Level 1 - may cause numerical issues
    Level1,
    /// Level 2 - usable for basic simulations
    Level2,
    /// Level 3 - suitable for most applications
    Level3,
    /// Level 4 - suitable for high-accuracy simulations
    Level4,
}

/// Statistical summary of a single metric
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MetricStatistics<T: RealField + Copy> {
    /// Minimum value
    pub min: T,
    /// Maximum value
    pub max: T,
    /// Mean value
    pub mean: T,
    /// Standard deviation
    pub std_dev: T,
}

/// Overall quality assessment
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityAssessment {
    /// Overall mesh quality
    pub overall_quality: QualityLevel,
    /// Aspect ratio quality
    pub aspect_ratio_quality: QualityLevel,
    /// Skewness quality
    pub skewness_quality: QualityLevel,
    /// Orthogonality quality
    pub orthogonality_quality: QualityLevel,
    /// Improvement recommendations
    pub recommendations: Vec<String>,
}

/// Trait for mesh quality assessment strategies
pub trait MeshQuality<T: RealField + Copy>: Send + Sync {
    /// Compute aspect ratio for a specific element
    fn aspect_ratio(&self, mesh: &Mesh<T>, element_idx: usize) -> Result<T>;

    /// Compute skewness for a specific element
    fn skewness(&self, mesh: &Mesh<T>, element_idx: usize) -> Result<T>;

    /// Generate a comprehensive quality report for the entire mesh
    fn quality_report(&self, mesh: &Mesh<T>) -> Result<QualityReport<T>>;
}

/// Comprehensive mesh quality report
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityReport<T: RealField + Copy> {
    /// Overall quality score (0.0 to 1.0, where 1.0 is perfect)
    pub overall_score: T,
    /// Detailed statistics for various metrics
    pub statistics: QualityStatistics<T>,
    /// Indices of elements that fail quality thresholds
    pub problematic_elements: Vec<usize>,
    /// Automated recommendations for improving mesh quality
    pub recommendations: Vec<String>,
}

/// Statistical summary of mesh quality metrics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityStatistics<T: RealField + Copy> {
    /// Minimum aspect ratio found
    pub min_aspect_ratio: T,
    /// Maximum aspect ratio found
    pub max_aspect_ratio: T,
    /// Average aspect ratio across all elements
    pub avg_aspect_ratio: T,
    /// Minimum skewness found
    pub min_skewness: T,
    /// Maximum skewness found
    pub max_skewness: T,
    /// Average skewness across all elements
    pub avg_skewness: T,
    /// Distribution of quality scores into buckets
    pub quality_distribution: HashMap<String, usize>,
}

#[cfg(test)]
mod tests {
    use super::{MeshQualityService, MetricStatistics, QualityLevel};

    fn stats(min: f64, max: f64) -> MetricStatistics<f64> {
        MetricStatistics {
            min,
            max,
            mean: (min + max) * 0.5,
            std_dev: 0.0,
        }
    }

    #[test]
    fn mesh_quality_reports_level4_when_all_metrics_are_inside_best_thresholds() {
        let assessment = MeshQualityService::assess_quality(
            &stats(1.0, 1.5),
            &stats(0.0, 0.2),
            &stats(0.97, 1.0),
        );

        assert_eq!(assessment.aspect_ratio_quality, QualityLevel::Level4);
        assert_eq!(assessment.skewness_quality, QualityLevel::Level4);
        assert_eq!(assessment.orthogonality_quality, QualityLevel::Level4);
        assert_eq!(assessment.overall_quality, QualityLevel::Level4);
        assert_eq!(
            assessment.recommendations,
            vec!["Mesh quality is acceptable for CFD simulation".to_string()]
        );
    }

    #[test]
    fn mesh_quality_reports_worst_metric_and_targeted_recommendations() {
        let assessment = MeshQualityService::assess_quality(
            &stats(1.0, 12.0),
            &stats(0.0, 0.9),
            &stats(0.6, 1.0),
        );

        assert_eq!(assessment.aspect_ratio_quality, QualityLevel::Level1);
        assert_eq!(assessment.skewness_quality, QualityLevel::Level1);
        assert_eq!(assessment.orthogonality_quality, QualityLevel::Level1);
        assert_eq!(assessment.overall_quality, QualityLevel::Level1);
        assert_eq!(
            assessment.recommendations,
            vec![
                "Consider refining mesh in regions with high aspect ratio cells".to_string(),
                "Improve mesh quality by reducing cell skewness".to_string(),
                "Enhance mesh orthogonality for better numerical accuracy".to_string(),
            ]
        );
    }

    #[test]
    fn mesh_quality_threshold_boundaries_use_strict_comparison_contracts() {
        let assessment = MeshQualityService::assess_quality(
            &stats(1.0, 5.0),
            &stats(0.0, 0.5),
            &stats(0.85, 1.0),
        );

        assert_eq!(assessment.aspect_ratio_quality, QualityLevel::Level2);
        assert_eq!(assessment.skewness_quality, QualityLevel::Level2);
        assert_eq!(assessment.orthogonality_quality, QualityLevel::Level2);
        assert_eq!(assessment.overall_quality, QualityLevel::Level2);
    }
}
