//! Geometric domain representations for 1D, 2D, and 3D space
//!
//! This module provides the fundamental building blocks for defining the
//! computational space in which simulations are performed.

use nalgebra::{Point1, Point2, Point3, RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Trait for geometric shapes and spatial operations
pub trait Geometry<T: RealField + Copy>: Send + Sync {
    /// Compute distance between two points
    fn distance(&self, p1: &Point3<T>, p2: &Point3<T>) -> T;

    /// Compute normal vector at a point
    fn normal(&self, point: &Point3<T>) -> Vector3<T>;

    /// Check if a point is inside the geometry
    fn contains(&self, point: &Point3<T>) -> bool;

    /// Project a point onto the nearest surface of the geometry
    fn project(&self, point: &Point3<T>) -> Point3<T>;

    /// Compute axis-aligned bounding box
    fn bounding_box(&self) -> (Point3<T>, Point3<T>);
}

/// Helper function to order two values
#[inline]
pub fn order<T: RealField>(v1: T, v2: T) -> (T, T) {
    if v1 <= v2 {
        (v1, v2)
    } else {
        (v2, v1)
    }
}

/// Trait for computational domains
pub trait Domain<T: RealField + Copy>: Send + Sync {
    /// Get the dimensionality of the domain (1, 2, or 3)
    fn dimension(&self) -> usize;

    /// Get the volume (or area/length) of the domain
    fn volume(&self) -> T;

    /// Check if a point is inside the domain (dimension-specific)
    fn contains_1d(&self, _point: &Point1<T>) -> Option<bool> {
        None
    }

    /// Check if a point is inside the domain (dimension-specific)
    fn contains_2d(&self, _point: &Point2<T>) -> Option<bool> {
        None
    }

    /// Check if a point is inside the domain (dimension-specific)
    fn contains_3d(&self, _point: &Point3<T>) -> Option<bool> {
        None
    }
}

/// 1D domain (line segment)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain1D<T: RealField + Copy> {
    /// Start point
    pub start: T,
    /// End point
    pub end: T,
}

impl<T: RealField + Copy> Domain1D<T> {
    /// Create a new 1D domain with automatic ordering
    pub fn new(p1: T, p2: T) -> Self {
        let (start, end) = order(p1, p2);
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

    /// Check if a point is within the domain
    pub fn contains(&self, point: &Point1<T>) -> bool {
        point.x >= self.start && point.x <= self.end
    }
}

impl<T: RealField + Copy> Domain<T> for Domain1D<T> {
    fn dimension(&self) -> usize {
        1
    }

    fn volume(&self) -> T {
        self.length()
    }

    fn contains_1d(&self, point: &Point1<T>) -> Option<bool> {
        Some(self.contains(point))
    }
}

/// 2D rectangular domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain2D<T: RealField + Copy> {
    /// Minimum corner
    pub min: Point2<T>,
    /// Maximum corner
    pub max: Point2<T>,
}

impl<T: RealField + Copy> Domain2D<T> {
    /// Create a new 2D domain with automatic ordering
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self {
        let (x_min, x_max) = order(p1.x, p2.x);
        let (y_min, y_max) = order(p1.y, p2.y);
        Self {
            min: Point2::new(x_min, y_min),
            max: Point2::new(x_max, y_max),
        }
    }

    /// Create from scalar coordinates
    pub fn from_scalars(x1: T, y1: T, x2: T, y2: T) -> Self {
        Self::new(Point2::new(x1, y1), Point2::new(x2, y2))
    }

    /// Get the center of the domain
    pub fn center(&self) -> Point2<T> {
        let two = T::one() + T::one();
        Point2::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two,
        )
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

    /// Check if a point is within the domain
    pub fn contains(&self, point: &Point2<T>) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
    }
}

impl<T: RealField + Copy> Domain<T> for Domain2D<T> {
    fn dimension(&self) -> usize {
        2
    }

    fn volume(&self) -> T {
        self.area()
    }

    fn contains_2d(&self, point: &Point2<T>) -> Option<bool> {
        Some(self.contains(point))
    }
}

/// 3D box domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain3D<T: RealField + Copy> {
    /// Minimum corner
    pub min: Point3<T>,
    /// Maximum corner
    pub max: Point3<T>,
}

impl<T: RealField + Copy> Domain3D<T> {
    /// Create a new 3D domain with automatic ordering
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self {
        let (x_min, x_max) = order(p1.x, p2.x);
        let (y_min, y_max) = order(p1.y, p2.y);
        let (z_min, z_max) = order(p1.z, p2.z);
        Self {
            min: Point3::new(x_min, y_min, z_min),
            max: Point3::new(x_max, y_max, z_max),
        }
    }

    /// Create from scalar coordinates
    pub fn from_scalars(x1: T, y1: T, z1: T, x2: T, y2: T, z2: T) -> Self {
        Self::new(Point3::new(x1, y1, z1), Point3::new(x2, y2, z2))
    }

    /// Get the center of the domain
    pub fn center(&self) -> Point3<T> {
        let two = T::one() + T::one();
        Point3::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two,
            (self.min.z + self.max.z) / two,
        )
    }

    /// Check if a point is within the domain
    pub fn contains(&self, point: &Point3<T>) -> bool {
        point.x >= self.min.x
            && point.x <= self.max.x
            && point.y >= self.min.y
            && point.y <= self.max.y
            && point.z >= self.min.z
            && point.z <= self.max.z
    }

    /// Get volume
    pub fn volume(&self) -> T {
        (self.max.x - self.min.x) * (self.max.y - self.min.y) * (self.max.z - self.min.z)
    }
}

impl<T: RealField + Copy> Domain<T> for Domain3D<T> {
    fn dimension(&self) -> usize {
        3
    }

    fn volume(&self) -> T {
        self.volume()
    }

    fn contains_3d(&self, point: &Point3<T>) -> Option<bool> {
        Some(self.contains(point))
    }
}

/// Generic domain that can be 1D, 2D, or 3D
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum AnyDomain<T: RealField + Copy> {
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

    fn volume(&self) -> T {
        match self {
            Self::D1(d) => d.volume(),
            Self::D2(d) => d.volume(),
            Self::D3(d) => d.volume(),
        }
    }

    fn contains_1d(&self, point: &Point1<T>) -> Option<bool> {
        match self {
            Self::D1(d) => d.contains_1d(point),
            _ => None,
        }
    }

    fn contains_2d(&self, point: &Point2<T>) -> Option<bool> {
        match self {
            Self::D2(d) => d.contains_2d(point),
            _ => None,
        }
    }

    fn contains_3d(&self, point: &Point3<T>) -> Option<bool> {
        match self {
            Self::D3(d) => d.contains_3d(point),
            _ => None,
        }
    }
}

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
