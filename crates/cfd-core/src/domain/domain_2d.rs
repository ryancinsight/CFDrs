//! 2D domain implementation

use super::common::{order, Domain};
use nalgebra::{Point2, RealField};
use serde::{Deserialize, Serialize};

/// 2D rectangular domain
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain2D<T: RealField + Copy> {
    /// Minimum corner
    pub min: Point2<T>,
    /// Maximum corner
    pub max: Point2<T>,
}

impl<T: RealField + Copy> Domain2D<T> {
    /// Create a new 2D domain from corner points.
    /// Coordinates are automatically ordered to ensure min <= max on each axis.
    pub fn new(p1: Point2<T>, p2: Point2<T>) -> Self {
        let (x_min, x_max) = order(p1.x, p2.x);
        let (y_min, y_max) = order(p1.y, p2.y);
        Self {
            min: Point2::new(x_min, y_min),
            max: Point2::new(x_max, y_max),
        }
    }

    /// Get the center of the domain
    pub fn center(&self) -> Point2<T> {
        let two = T::one() + T::one();
        Point2::new(
            (self.min.x + self.max.x) / two,
            (self.min.y + self.max.y) / two,
        )
    }

    /// Create a new 2D domain from scalar coordinates.
    /// Coordinates are automatically ordered to ensure min <= max on each axis.
    pub fn from_scalars(x1: T, y1: T, x2: T, y2: T) -> Self {
        Self::new(Point2::new(x1, y1), Point2::new(x2, y2))
    }

    /// Create from two points (alias for new)
    pub fn from_points(p1: Point2<T>, p2: Point2<T>) -> Self {
        Self::new(p1, p2)
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
        Some(
            point.x >= self.min.x
                && point.x <= self.max.x
                && point.y >= self.min.y
                && point.y <= self.max.y,
        )
    }
}
