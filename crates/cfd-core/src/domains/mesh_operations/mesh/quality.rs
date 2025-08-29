//! Mesh quality metrics and analysis

use super::types::{Element, Mesh};
use crate::domains::mesh_operations::element::ElementType;
use crate::error::{Error, Result};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};

/// Quality metrics for mesh analysis
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityMetrics<T: RealField + Copy> {
    /// Minimum element quality
    pub min_quality: T,
    /// Maximum element quality
    pub max_quality: T,
    /// Average element quality
    pub avg_quality: T,
    /// Aspect ratio statistics
    pub aspect_ratio: AspectRatioStats<T>,
    /// Skewness statistics
    pub skewness: SkewnessStats<T>,
}

/// Aspect ratio statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AspectRatioStats<T: RealField + Copy> {
    pub min: T,
    pub max: T,
    pub avg: T,
}

/// Skewness statistics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SkewnessStats<T: RealField + Copy> {
    pub min: T,
    pub max: T,
    pub avg: T,
}

/// Trait for mesh quality analysis
pub trait MeshQuality<T: RealField + Copy> {
    /// Compute quality metrics
    fn compute_quality(&self) -> Result<QualityMetrics<T>>;

    /// Check element quality
    fn check_element_quality(&self, element: &Element) -> T;

    /// Compute element volume
    fn element_volume(&self, element: &Element) -> T;

    /// Compute element aspect ratio
    fn element_aspect_ratio(&self, element: &Element) -> T;
}

impl<T: RealField + Copy + Float + FromPrimitive> MeshQuality<T> for Mesh<T> {
    fn compute_quality(&self) -> Result<QualityMetrics<T>> {
        if self.elements.is_empty() {
            return Err(Error::InvalidInput("Mesh has no elements".into()));
        }

        let mut qualities = Vec::with_capacity(self.elements.len());
        let mut aspect_ratios = Vec::with_capacity(self.elements.len());

        for element in &self.elements {
            qualities.push(self.check_element_quality(element));
            aspect_ratios.push(self.element_aspect_ratio(element));
        }

        let sum_quality: T = qualities.iter().copied().fold(T::zero(), |a, b| a + b);
        let sum_aspect: T = aspect_ratios.iter().copied().fold(T::zero(), |a, b| a + b);
        let count = T::from_usize(qualities.len()).unwrap_or_else(T::one);

        Ok(QualityMetrics {
            min_quality: *qualities
                .iter()
                .min_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap(),
            max_quality: *qualities
                .iter()
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .unwrap(),
            avg_quality: sum_quality / count,
            aspect_ratio: AspectRatioStats {
                min: *aspect_ratios
                    .iter()
                    .min_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap(),
                max: *aspect_ratios
                    .iter()
                    .max_by(|a, b| a.partial_cmp(b).unwrap())
                    .unwrap(),
                avg: sum_aspect / count,
            },
            skewness: SkewnessStats {
                min: T::zero(),
                max: T::one(),
                avg: T::from_f64(0.5).unwrap_or_else(|| {
                    T::one() / T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one())
                }),
            },
        })
    }

    fn check_element_quality(&self, element: &Element) -> T {
        // Simplified quality metric based on element type
        match element.element_type {
            ElementType::Triangle => self.triangle_quality(element),
            ElementType::Tetrahedron => self.tetrahedron_quality(element),
            _ => T::one(),
        }
    }

    fn element_volume(&self, element: &Element) -> T {
        match element.element_type {
            ElementType::Tetrahedron if element.nodes.len() == 4 => {
                let p0 = self.nodes[element.nodes[0]];
                let p1 = self.nodes[element.nodes[1]];
                let p2 = self.nodes[element.nodes[2]];
                let p3 = self.nodes[element.nodes[3]];

                let v1 = p1 - p0;
                let v2 = p2 - p0;
                let v3 = p3 - p0;

                Float::abs(v1.cross(&v2).dot(&v3) / T::from_f64(6.0).unwrap_or_else(T::one))
            }
            _ => T::zero(),
        }
    }

    fn element_aspect_ratio(&self, _element: &Element) -> T {
        // Simplified aspect ratio calculation
        T::one()
    }
}

impl<T: RealField + Copy + Float + FromPrimitive> Mesh<T> {
    fn triangle_quality(&self, element: &Element) -> T {
        if element.nodes.len() != 3 {
            return T::zero();
        }

        let p0 = self.nodes[element.nodes[0]];
        let p1 = self.nodes[element.nodes[1]];
        let p2 = self.nodes[element.nodes[2]];

        let a = (p1 - p0).norm();
        let b = (p2 - p1).norm();
        let c = (p0 - p2).norm();

        let two = T::from_f64(2.0).unwrap_or_else(|| T::one() + T::one());
        let four = T::from_f64(4.0).unwrap_or_else(|| two + two);
        let s = (a + b + c) / two;
        let area = Float::sqrt(s * (s - a) * (s - b) * (s - c));

        if area > T::zero() {
            let radius_ratio = area / (a * b * c / four);
            Float::min(radius_ratio, T::one())
        } else {
            T::zero()
        }
    }

    fn tetrahedron_quality(&self, element: &Element) -> T {
        if element.nodes.len() != 4 {
            return T::zero();
        }

        let volume = self.element_volume(element);
        if volume > T::zero() {
            T::one() // Simplified for now
        } else {
            T::zero()
        }
    }
}
