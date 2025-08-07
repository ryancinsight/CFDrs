//! Computational domain representations.

use nalgebra::{Point3, RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Trait for computational domains
pub trait Domain<T: RealField>: Send + Sync {
    /// Get the dimensionality of the domain (1, 2, or 3)
    fn dimension(&self) -> usize;

    /// Check if a point is inside the domain
    fn contains(&self, point: &Point3<T>) -> bool;

    /// Get the bounding box of the domain
    fn bounding_box(&self) -> (Point3<T>, Point3<T>);

    /// Get the volume (or area/length) of the domain
    fn volume(&self) -> T;
}

/// 1D domain (line segment)
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Domain1D<T: RealField> {
    /// Start point
    pub start: T,
    /// End point
    pub end: T,
}

impl<T: RealField> Domain1D<T> {
    /// Create a new 1D domain
    pub const fn new(start: T, end: T) -> Self {
        Self { start, end }
    }

    /// Get the length of the domain
    pub fn length(&self) -> T {
        (self.end.clone() - self.start.clone()).abs()
    }

    /// Get the center of the domain
    pub fn center(&self) -> T {
        let two = T::one() + T::one();
        (self.start.clone() + self.end.clone()) / two
    }
}

impl<T: RealField> Domain<T> for Domain1D<T> {
    fn dimension(&self) -> usize {
        1
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        let x = point.x.clone();
        let min = self.start.clone().min(self.end.clone());
        let max = self.start.clone().max(self.end.clone());
        x >= min && x <= max
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        let min = self.start.clone().min(self.end.clone());
        let max = self.start.clone().max(self.end.clone());
        (
            Point3::new(min, T::zero(), T::zero()),
            Point3::new(max, T::zero(), T::zero()),
        )
    }

    fn volume(&self) -> T {
        self.length()
    }
}

/// 2D rectangular domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain2D<T: RealField> {
    /// Minimum corner
    pub min: Point3<T>,
    /// Maximum corner
    pub max: Point3<T>,
}

impl<T: RealField> Domain2D<T> {
    /// Create a new 2D domain
    pub fn new(x_min: T, y_min: T, x_max: T, y_max: T) -> Self {
        Self {
            min: Point3::new(x_min, y_min, T::zero()),
            max: Point3::new(x_max, y_max, T::zero()),
        }
    }

    /// Get the width of the domain
    pub fn width(&self) -> T {
        self.max.x.clone() - self.min.x.clone()
    }

    /// Get the height of the domain
    pub fn height(&self) -> T {
        self.max.y.clone() - self.min.y.clone()
    }

    /// Get the area of the domain
    pub fn area(&self) -> T {
        self.width() * self.height()
    }
}

impl<T: RealField> Domain<T> for Domain2D<T> {
    fn dimension(&self) -> usize {
        2
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        (self.min.clone(), self.max.clone())
    }

    fn volume(&self) -> T {
        self.area()
    }
}

/// 3D box domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain3D<T: RealField> {
    /// Minimum corner
    pub min: Point3<T>,
    /// Maximum corner
    pub max: Point3<T>,
}

impl<T: RealField> Domain3D<T> {
    /// Create a new 3D domain
    pub const fn new(min: Point3<T>, max: Point3<T>) -> Self {
        Self { min, max }
    }

    /// Get the width (x dimension)
    pub fn width(&self) -> T {
        self.max.x.clone() - self.min.x.clone()
    }

    /// Get the height (y dimension)
    pub fn height(&self) -> T {
        self.max.y.clone() - self.min.y.clone()
    }

    /// Get the depth (z dimension)
    pub fn depth(&self) -> T {
        self.max.z.clone() - self.min.z.clone()
    }

    /// Get the center of the domain
    pub fn center(&self) -> Point3<T> {
        let two = T::one() + T::one();
        Point3::new(
            (self.min.x.clone() + self.max.x.clone()) / two.clone(),
            (self.min.y.clone() + self.max.y.clone()) / two.clone(),
            (self.min.z.clone() + self.max.z.clone()) / two,
        )
    }

    /// Create from center and half-extents
    pub fn from_center_half_extents(center: Point3<T>, half_extents: Vector3<T>) -> Self {
        Self {
            min: center.clone() - half_extents.clone(),
            max: center + half_extents,
        }
    }

    /// Get the diagonal vector
    pub fn diagonal(&self) -> Vector3<T> {
        self.max.clone() - self.min.clone()
    }

    /// Get the volume
    pub fn volume(&self) -> T {
        let dims = self.diagonal();
        dims.x.clone() * dims.y.clone() * dims.z.clone()
    }
}

impl<T: RealField> Domain<T> for Domain3D<T> {
    fn dimension(&self) -> usize {
        3
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        (self.min.clone(), self.max.clone())
    }

    fn volume(&self) -> T {
        let dx = self.max.x.clone() - self.min.x.clone();
        let dy = self.max.y.clone() - self.min.y.clone();
        let dz = self.max.z.clone() - self.min.z.clone();
        dx * dy * dz
    }
}

/// Generic domain that can be 1D, 2D, or 3D
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AnyDomain<T: RealField> {
    /// 1D domain
    D1(Domain1D<T>),
    /// 2D domain
    D2(Domain2D<T>),
    /// 3D domain
    D3(Domain3D<T>),
}

impl<T: RealField> Domain<T> for AnyDomain<T> {
    fn dimension(&self) -> usize {
        match self {
            Self::D1(d) => d.dimension(),
            Self::D2(d) => d.dimension(),
            Self::D3(d) => d.dimension(),
        }
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        match self {
            Self::D1(d) => d.contains(point),
            Self::D2(d) => d.contains(point),
            Self::D3(d) => d.contains(point),
        }
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        match self {
            Self::D1(d) => d.bounding_box(),
            Self::D2(d) => d.bounding_box(),
            Self::D3(d) => d.bounding_box(),
        }
    }

    fn volume(&self) -> T {
        match self {
            Self::D1(d) => d.volume(),
            Self::D2(d) => d.volume(),
            Self::D3(d) => d.volume(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_domain_1d() {
        let domain = Domain1D::new(0.0, 1.0);
        assert_eq!(domain.dimension(), 1);
        assert_relative_eq!(domain.length(), 1.0);
        assert!(domain.contains(&Point3::new(0.5, 0.0, 0.0)));
        assert!(!domain.contains(&Point3::new(1.5, 0.0, 0.0)));
    }

    #[test]
    fn test_domain_2d() {
        let domain = Domain2D::new(0.0, 0.0, 2.0, 3.0);
        assert_eq!(domain.dimension(), 2);
        assert_relative_eq!(domain.area(), 6.0);
        assert!(domain.contains(&Point3::new(1.0, 1.0, 0.0)));
        assert!(!domain.contains(&Point3::new(3.0, 1.0, 0.0)));
    }

    #[test]
    fn test_domain_3d() {
        let domain = Domain3D::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 4.0),
        );
        assert_eq!(domain.dimension(), 3);
        assert_relative_eq!(domain.volume(), 24.0);
        assert!(domain.contains(&Point3::new(1.0, 1.0, 1.0)));
        assert!(!domain.contains(&Point3::new(3.0, 1.0, 1.0)));
    }
}