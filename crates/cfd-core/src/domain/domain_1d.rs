//! 1D domain implementation

use nalgebra::{Point1, RealField};
use serde::{Deserialize, Serialize};
use super::common::Domain;

/// 1D domain (line segment)
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Domain1D<T: RealField + Copy> {
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
        Some(point.x >= self.start && point.x <= self.end)
    }
}