//! Common domain traits and utilities

use nalgebra::{Point1, Point2, Point3, RealField};

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
    /// Returns Some(true) if contained, Some(false) if not contained,
    /// None if called on an incompatible dimension
    fn contains_1d(&self, _point: &Point1<T>) -> Option<bool> {
        None // Default: not supported for this dimension
    }

    /// Check if a point is inside the domain (dimension-specific)
    /// Returns Some(true) if contained, Some(false) if not contained,
    /// None if called on an incompatible dimension
    fn contains_2d(&self, _point: &Point2<T>) -> Option<bool> {
        None // Default: not supported for this dimension
    }

    /// Check if a point is inside the domain (dimension-specific)
    /// Returns Some(true) if contained, Some(false) if not contained,
    /// None if called on an incompatible dimension
    fn contains_3d(&self, _point: &Point3<T>) -> Option<bool> {
        None // Default: not supported for this dimension
    }
}