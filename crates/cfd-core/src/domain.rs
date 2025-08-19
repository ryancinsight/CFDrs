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
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain1D<T: RealField> {
    /// Start point (guaranteed to be <= end after construction)
    pub start: T,
    /// End point (guaranteed to be >= start after construction)
    pub end: T,
}

impl<T: RealField + Copy> Domain1D<T> {
    /// Create a new 1D domain. The start and end points are automatically ordered
    /// so that start <= end.
    pub fn new(p1: T, p2: T) -> Self {
        let (start, end) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
        Self { start, end }
    }

    /// Get the length of the domain
    pub fn length(&self) -> T {
        self.end - self.start
    }

    /// Get the center of the domain
    pub fn center(&self) -> T {
        let two = T::one() + T::one();
        (self.start + self.end) / two
    }
}

impl<T: RealField + Copy> Domain<T> for Domain1D<T> {
    fn dimension(&self) -> usize {
        1
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        // Since we enforce start <= end in the constructor, we can simplify this
        point.x >= self.start && point.x <= self.end
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        // Since we enforce start <= end in the constructor, we can simplify this
        (
            Point3::new(self.start, T::zero(), T::zero()),
            Point3::new(self.end, T::zero(), T::zero()),
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

impl<T: RealField + Copy> Domain2D<T> {
    /// Create a new 2D domain from scalar coordinates.
    pub fn new(x_min: T, y_min: T, x_max: T, y_max: T) -> Self {
        Self {
            min: Point3::new(x_min, y_min, T::zero()),
            max: Point3::new(x_max, y_max, T::zero()),
        }
    }

    /// Create a new 2D domain from corner points.
    pub fn from_points(min: Point3<T>, max: Point3<T>) -> Self {
        // In debug builds, ensure z-components are zero for 2D domains
        debug_assert!(min.z.is_zero());
        debug_assert!(max.z.is_zero());
        Self { min, max }
    }

    /// Get the width of the domain
    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    /// Get the height of the domain
    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    /// Get the area of the domain
    pub fn area(&self) -> T {
        self.width() * self.height()
    }
}

impl<T: RealField + Copy> Domain<T> for Domain2D<T> {
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

impl<T: RealField + Copy> Domain3D<T> {
    /// Create a new 3D domain
    pub const fn new(min: Point3<T>, max: Point3<T>) -> Self {
        Self { min, max }
    }

    /// Get the width (x dimension)
    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    /// Get the height (y dimension)
    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    /// Get the depth (z dimension)
    pub fn depth(&self) -> T {
        self.max.z - self.min.z
    }

    /// Get the center of the domain
    pub fn center(&self) -> Point3<T> {
        let two = T::one() + T::one();
        Point3::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two.clone(),
            (self.min.z + self.max.z) / two.clone(),
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
        dims.x * dims.y * dims.z
    }
}

impl<T: RealField + Copy> Domain<T> for Domain3D<T> {
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
        // Delegate to the inherent volume method to avoid code duplication
        Domain3D::volume(self)
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

impl<T: RealField + Copy> Domain<T> for AnyDomain<T> {
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

// Ergonomic From implementations for AnyDomain conversions
impl<T: RealField + Copy> From<Domain1D<T>> for AnyDomain<T> {
    fn from(domain: Domain1D<T>) -> Self {
        AnyDomain::D1(domain)
    }
}

impl<T: RealField + Copy> From<Domain2D<T>> for AnyDomain<T> {
    fn from(domain: Domain2D<T>) -> Self {
        AnyDomain::D2(domain)
    }
}

impl<T: RealField + Copy> From<Domain3D<T>> for AnyDomain<T> {
    fn from(domain: Domain3D<T>) -> Self {
        AnyDomain::D3(domain)
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
        
        // Test automatic ordering
        let domain_reversed = Domain1D::new(1.0, 0.0);
        assert_eq!(domain_reversed.start, 0.0);
        assert_eq!(domain_reversed.end, 1.0);
        assert_relative_eq!(domain_reversed.length(), 1.0);
    }

    #[test]
    fn test_domain_2d() {
        let domain = Domain2D::new(0.0, 0.0, 2.0, 3.0);
        assert_eq!(domain.dimension(), 2);
        assert_relative_eq!(domain.area(), 6.0);
        assert!(domain.contains(&Point3::new(1.0, 1.0, 0.0)));
        assert!(!domain.contains(&Point3::new(3.0, 1.0, 0.0)));

        // Test from_points constructor
        let domain2 = Domain2D::from_points(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 3.0, 0.0),
        );
        assert_relative_eq!(domain2.area(), 6.0);
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

    #[test]
    fn test_any_domain_from_conversions() {
        let domain_1d = Domain1D::new(0.0, 1.0);
        let any_domain: AnyDomain<f64> = domain_1d.into();
        assert_eq!(any_domain.dimension(), 1);

        let domain_2d = Domain2D::new(0.0, 0.0, 1.0, 1.0);
        let any_domain: AnyDomain<f64> = domain_2d.into();
        assert_eq!(any_domain.dimension(), 2);

        let domain_3d = Domain3D::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        );
        let any_domain: AnyDomain<f64> = domain_3d.into();
        assert_eq!(any_domain.dimension(), 3);
    }
}