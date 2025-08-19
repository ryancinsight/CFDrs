//! Fixed computational domain representations with proper type safety
//!
//! This module provides corrected domain implementations addressing:
//! - Type safety for 2D domains using Point2
//! - Consistent const fn usage
//! - Proper validation of invariants
//! - Ergonomic enum accessors

use nalgebra::{Point2, Point3, RealField, Vector2, Vector3};
use serde::{Deserialize, Serialize};

/// Trait for computational domains with proper const bounds
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

/// 1D domain (line segment) with validation
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain1D<T: RealField> {
    /// Start point (guaranteed to be <= end after construction)
    start: T,
    /// End point (guaranteed to be >= start after construction)
    end: T,
}

impl<T: RealField> Domain1D<T> {
    /// Create a new 1D domain with automatic ordering
    /// 
    /// # Invariant
    /// Ensures start <= end by swapping if necessary
    pub fn new(p1: T, p2: T) -> Self {
        let (start, end) = if p1 <= p2 { (p1, p2) } else { (p2, p1) };
        Self { start, end }
    }

    /// Create a validated 1D domain, returning error if invalid
    pub fn new_validated(start: T, end: T) -> Result<Self, DomainError> {
        if start > end {
            return Err(DomainError::InvalidBounds {
                message: "Start must be <= end".into(),
            });
        }
        Ok(Self { start, end })
    }

    /// Get the length of the domain
    #[inline]
    pub fn length(&self) -> T {
        self.end - self.start
    }

    /// Get the center of the domain
    #[inline]
    pub fn center(&self) -> T {
        let two = T::one() + T::one();
        (self.start + self.end) / two
    }
}

/// 2D rectangular domain with proper type safety
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain2D<T: RealField> {
    /// Minimum corner in 2D space
    min: Point2<T>,
    /// Maximum corner in 2D space
    max: Point2<T>,
}

impl<T: RealField> Domain2D<T> {
    /// Create a new 2D domain from scalar coordinates
    /// 
    /// Now const fn for consistency with Domain3D!
    pub const fn new_from_scalars(x_min: T, y_min: T, x_max: T, y_max: T) -> Self {
        Self {
            min: Point2::new(x_min, y_min),
            max: Point2::new(x_max, y_max),
        }
    }

    /// Create a new 2D domain from points
    pub const fn new(min: Point2<T>, max: Point2<T>) -> Self {
        Self { min, max }
    }

    /// Create a validated 2D domain
    pub fn new_validated(min: Point2<T>, max: Point2<T>) -> Result<Self, DomainError> {
        if min.x > max.x || min.y > max.y {
            return Err(DomainError::InvalidBounds {
                message: "Min must be <= max in all dimensions".into(),
            });
        }
        Ok(Self { min, max })
    }

    /// Get the width of the domain
    #[inline]
    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    /// Get the height of the domain
    #[inline]
    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    /// Get the area of the domain
    #[inline]
    pub fn area(&self) -> T {
        self.width() * self.height()
    }

    /// Get the center of the domain
    #[inline]
    pub fn center(&self) -> Point2<T> {
        let two = T::one() + T::one();
        Point2::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two,
        )
    }
}

impl<T: RealField> Domain<T> for Domain2D<T> {
    fn dimension(&self) -> usize {
        2
    }

    fn contains(&self, point: &Point3<T>) -> bool {
        // Properly check z-component for 2D domain
        if !point.z.is_zero() {
            return false;
        }
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }

    fn bounding_box(&self) -> (Point3<T>, Point3<T>) {
        // Lift 2D points to 3D space with z=0
        (
            Point3::new(self.min.x, self.min.y, T::zero()),
            Point3::new(self.max.x, self.max.y, T::zero()),
        )
    }

    fn volume(&self) -> T {
        self.area()
    }
}

/// 3D box domain with validation
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain3D<T: RealField> {
    /// Minimum corner
    min: Point3<T>,
    /// Maximum corner
    max: Point3<T>,
}

impl<T: RealField> Domain3D<T> {
    /// Create a new 3D domain (const fn for performance)
    pub const fn new(min: Point3<T>, max: Point3<T>) -> Self {
        Self { min, max }
    }

    /// Create a validated 3D domain
    pub fn new_validated(min: Point3<T>, max: Point3<T>) -> Result<Self, DomainError> {
        if min.x > max.x || min.y > max.y || min.z > max.z {
            return Err(DomainError::InvalidBounds {
                message: "Min must be <= max in all dimensions".into(),
            });
        }
        Ok(Self { min, max })
    }

    /// Get the width (x dimension)
    #[inline]
    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    /// Get the height (y dimension)
    #[inline]
    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    /// Get the depth (z dimension)
    #[inline]
    pub fn depth(&self) -> T {
        self.max.z - self.min.z
    }

    /// Get the diagonal vector
    #[inline]
    pub fn diagonal(&self) -> Vector3<T> {
        self.max - self.min
    }

    /// Get the center of the domain
    #[inline]
    pub fn center(&self) -> Point3<T> {
        let two = T::one() + T::one();
        Point3::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two,
            (self.min.z + self.max.z) / two,
        )
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
        (self.min, self.max)
    }

    fn volume(&self) -> T {
        // Direct implementation without delegation
        let dims = self.max - self.min;
        dims.x * dims.y * dims.z
    }
}

/// Generic domain enum with ergonomic accessors
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AnyDomain<T: RealField> {
    /// 1D domain
    D1(Domain1D<T>),
    /// 2D domain
    D2(Domain2D<T>),
    /// 3D domain
    D3(Domain3D<T>),
}

impl<T: RealField> AnyDomain<T> {
    /// Returns this domain as a `Domain1D` if it is one
    pub fn as_1d(&self) -> Option<&Domain1D<T>> {
        match self {
            Self::D1(d) => Some(d),
            _ => None,
        }
    }

    /// Returns this domain as a mutable `Domain1D` if it is one
    pub fn as_1d_mut(&mut self) -> Option<&mut Domain1D<T>> {
        match self {
            Self::D1(d) => Some(d),
            _ => None,
        }
    }

    /// Converts into `Domain1D` if it is one
    pub fn into_1d(self) -> Option<Domain1D<T>> {
        match self {
            Self::D1(d) => Some(d),
            _ => None,
        }
    }

    /// Returns this domain as a `Domain2D` if it is one
    pub fn as_2d(&self) -> Option<&Domain2D<T>> {
        match self {
            Self::D2(d) => Some(d),
            _ => None,
        }
    }

    /// Returns this domain as a mutable `Domain2D` if it is one
    pub fn as_2d_mut(&mut self) -> Option<&mut Domain2D<T>> {
        match self {
            Self::D2(d) => Some(d),
            _ => None,
        }
    }

    /// Converts into `Domain2D` if it is one
    pub fn into_2d(self) -> Option<Domain2D<T>> {
        match self {
            Self::D2(d) => Some(d),
            _ => None,
        }
    }

    /// Returns this domain as a `Domain3D` if it is one
    pub fn as_3d(&self) -> Option<&Domain3D<T>> {
        match self {
            Self::D3(d) => Some(d),
            _ => None,
        }
    }

    /// Returns this domain as a mutable `Domain3D` if it is one
    pub fn as_3d_mut(&mut self) -> Option<&mut Domain3D<T>> {
        match self {
            Self::D3(d) => Some(d),
            _ => None,
        }
    }

    /// Converts into `Domain3D` if it is one
    pub fn into_3d(self) -> Option<Domain3D<T>> {
        match self {
            Self::D3(d) => Some(d),
            _ => None,
        }
    }

    /// Get the dimensionality without matching
    pub fn dimensionality(&self) -> usize {
        match self {
            Self::D1(_) => 1,
            Self::D2(_) => 2,
            Self::D3(_) => 3,
        }
    }
}

impl<T: RealField> Domain<T> for AnyDomain<T> {
    fn dimension(&self) -> usize {
        self.dimensionality()
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

/// Domain-related errors
#[derive(Debug, Clone)]
pub enum DomainError {
    /// Invalid domain bounds
    InvalidBounds { message: String },
    /// Dimension mismatch
    DimensionMismatch { expected: usize, actual: usize },
}

impl std::fmt::Display for DomainError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidBounds { message } => write!(f, "Invalid domain bounds: {}", message),
            Self::DimensionMismatch { expected, actual } => {
                write!(f, "Dimension mismatch: expected {}, got {}", expected, actual)
            }
        }
    }
}

impl std::error::Error for DomainError {}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_domain_2d_type_safety() {
        let domain = Domain2D::new(
            Point2::new(0.0, 0.0),
            Point2::new(1.0, 1.0),
        );
        
        // This point has non-zero z, should not be contained
        let point_3d = Point3::new(0.5, 0.5, 0.1);
        assert!(!domain.contains(&point_3d));
        
        // This point has zero z, should be contained
        let point_2d = Point3::new(0.5, 0.5, 0.0);
        assert!(domain.contains(&point_2d));
    }

    #[test]
    fn test_domain_validation() {
        // Should fail validation
        let invalid = Domain3D::new_validated(
            Point3::new(10.0, 10.0, 10.0),
            Point3::new(0.0, 0.0, 0.0),
        );
        assert!(invalid.is_err());
        
        // Should pass validation
        let valid = Domain3D::new_validated(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(10.0, 10.0, 10.0),
        );
        assert!(valid.is_ok());
    }

    #[test]
    fn test_any_domain_accessors() {
        let domain_3d = Domain3D::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        );
        let any_domain = AnyDomain::D3(domain_3d);
        
        // Ergonomic access to specific type
        if let Some(d3) = any_domain.as_3d() {
            assert_relative_eq!(d3.depth(), 1.0);
        } else {
            panic!("Should be a 3D domain");
        }
        
        // Should return None for wrong type
        assert!(any_domain.as_2d().is_none());
    }

    #[test]
    fn test_const_fn_usage() {
        // This can now be used in const context!
        const DOMAIN_2D: Domain2D<f64> = Domain2D::new_from_scalars(0.0, 0.0, 1.0, 1.0);
        const DOMAIN_3D: Domain3D<f64> = Domain3D::new(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 1.0),
        );
        
        assert_eq!(DOMAIN_2D.dimension(), 2);
        assert_eq!(DOMAIN_3D.dimension(), 3);
    }
}