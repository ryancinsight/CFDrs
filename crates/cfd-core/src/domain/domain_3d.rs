//! 3D domain implementation

use nalgebra::{Point3, RealField, Vector3};
use serde::{Deserialize, Serialize};
use super::common::{Domain, order};

/// 3D box domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain3D<T: RealField + Copy> {
    /// Minimum corner
    pub min: Point3<T>,
    /// Maximum corner
    pub max: Point3<T>,
}

impl<T: RealField + Copy> Domain3D<T> {
    /// Create a new 3D domain from corner points.
    /// Coordinates are automatically ordered to ensure min <= max on each axis.
    pub fn new(p1: Point3<T>, p2: Point3<T>) -> Self {
        let (x_min, x_max) = order(p1.x, p2.x);
        let (y_min, y_max) = order(p1.y, p2.y);
        let (z_min, z_max) = order(p1.z, p2.z);
        Self {
            min: Point3::new(x_min, y_min, z_min),
            max: Point3::new(x_max, y_max, z_max),
        }
    }

    /// Create a new 3D domain from scalar coordinates.
    /// Coordinates are automatically ordered to ensure min <= max on each axis.
    pub fn from_scalars(x1: T, y1: T, z1: T, x2: T, y2: T, z2: T) -> Self {
        Self::new(Point3::new(x1, y1, z1), Point3::new(x2, y2, z2))
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
            (self.min.y + self.max.y) / two,
            (self.min.z + self.max.z) / two,
        )
    }

    /// Create from center and half-extents
    pub fn from_center_half_extents(center: Point3<T>, half_extents: Vector3<T>) -> Self {
        Self {
            min: center - half_extents,
            max: center + half_extents,
        }
    }

    /// Get the diagonal vector
    pub fn diagonal(&self) -> Vector3<T> {
        self.max - self.min
    }

    /// Get the volume
    pub fn volume(&self) -> T {
        self.width() * self.height() * self.depth()
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
}

impl<T: RealField + Copy> Domain<T> for Domain3D<T> {
    fn dimension(&self) -> usize {
        3
    }

    fn volume(&self) -> T {
        self.width() * self.height() * self.depth()
    }

    fn contains_3d(&self, point: &Point3<T>) -> Option<bool> {
        Some(
            point.x >= self.min.x
                && point.x <= self.max.x
                && point.y >= self.min.y
                && point.y <= self.max.y
                && point.z >= self.min.z
                && point.z <= self.max.z,
        )
    }
}