//! Vertex representation in 3D space

use nalgebra::{Point3, RealField};
use serde::{Deserialize, Serialize};

/// Vertex in 3D space
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Vertex<T: RealField + Copy> {
    /// Position in 3D space
    pub position: Point3<T>,
    /// Global ID for distributed meshes
    pub global_id: Option<usize>,
    /// Partition ID (rank) that owns this vertex
    pub partition_id: Option<usize>,
}

impl<T: RealField + Copy> Vertex<T> {
    /// Create a new vertex at the given position
    pub fn new(position: Point3<T>) -> Self {
        Self {
            position,
            global_id: None,
            partition_id: None,
        }
    }

    /// Create a new vertex from coordinates
    pub fn from_coords(x: T, y: T, z: T) -> Self {
        Self {
            position: Point3::new(x, y, z),
            global_id: None,
            partition_id: None,
        }
    }

    /// Set distributed mesh properties
    pub fn with_distributed_info(mut self, global_id: usize, partition_id: usize) -> Self {
        self.global_id = Some(global_id);
        self.partition_id = Some(partition_id);
        self
    }

    /// Distance to another vertex
    pub fn distance_to(&self, other: &Self) -> T {
        (self.position - other.position).norm()
    }
}
