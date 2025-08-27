//! Geometric entities for FVM

use nalgebra::{RealField, Vector2};

/// Face between two control volumes
#[derive(Debug, Clone)]
pub struct Face<T: RealField + Copy> {
    /// Face center position
    pub center: Vector2<T>,
    /// Face normal vector (unit)
    pub normal: Vector2<T>,
    /// Face area/length in 2D
    pub area: T,
    /// Owner cell index
    pub owner: usize,
    /// Neighbor cell index (None for boundary faces)
    pub neighbor: Option<usize>,
}

impl<T: RealField + Copy> Face<T> {
    /// Create a new face
    pub fn new(
        center: Vector2<T>,
        normal: Vector2<T>,
        area: T,
        owner: usize,
        neighbor: Option<usize>,
    ) -> Self {
        Self {
            center,
            normal: normal.normalize(),
            area,
            owner,
            neighbor,
        }
    }

    /// Check if this is a boundary face
    pub fn is_boundary(&self) -> bool {
        self.neighbor.is_none()
    }

    /// Get the flux through this face
    pub fn flux(&self, velocity: Vector2<T>) -> T {
        self.area * velocity.dot(&self.normal)
    }
}
