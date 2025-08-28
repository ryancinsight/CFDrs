//! Geometry trait for mesh operations

use nalgebra::{Point3, RealField, Vector3};

/// Geometry trait for mesh operations
pub trait Geometry<T: RealField + Copy>: Send + Sync {
    /// Compute distance between two points
    fn distance(&self, p1: &Point3<T>, p2: &Point3<T>) -> T;

    /// Compute normal vector at a point
    fn normal(&self, point: &Point3<T>) -> Vector3<T>;

    /// Check if a point is inside the geometry
    fn contains(&self, point: &Point3<T>) -> bool;

    /// Project a point onto the geometry
    fn project(&self, point: &Point3<T>) -> Point3<T>;

    /// Compute bounding box
    fn bounding_box(&self) -> (Point3<T>, Point3<T>);
}
