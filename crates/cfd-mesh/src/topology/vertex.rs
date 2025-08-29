//! Vertex representation in 3D space

use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};

/// Vertex in 3D space
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vertex<T: RealField + Copy> {
    /// Position in 3D space
    pub position: Point3<T>,
}

impl<T: RealField + Copy> Vertex<T> {
    /// Create a new vertex at the given position
    pub fn new(position: Point3<T>) -> Self {
        Self { position }
    }

    /// Create a new vertex from coordinates
    pub fn from_coords(x: T, y: T, z: T) -> Self {
        Self {
            position: Point3::new(x, y, z),
        }
    }

    /// Distance to another vertex
    pub fn distance_to(&self, other: &Self) -> T {
        (self.position - other.position).norm()
    }
}
