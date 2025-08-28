//! Lagrangian point representation for IBM

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Lagrangian point representing an immersed boundary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LagrangianPoint<T: RealField + Copy> {
    /// Position of the Lagrangian point
    pub position: Vector3<T>,
    /// Velocity at the Lagrangian point
    pub velocity: Vector3<T>,
    /// Force at the Lagrangian point
    pub force: Vector3<T>,
    /// Area/volume associated with this point
    pub weight: T,
}

impl<T: RealField + Copy> LagrangianPoint<T> {
    /// Create a new Lagrangian point
    pub fn new(position: Vector3<T>, weight: T) -> Self {
        Self {
            position,
            velocity: Vector3::zeros(),
            force: Vector3::zeros(),
            weight,
        }
    }

    /// Update position based on velocity
    pub fn update_position(&mut self, dt: T) {
        self.position += self.velocity * dt;
    }

    /// Reset forces to zero
    pub fn reset_force(&mut self) {
        self.force = Vector3::zeros();
    }
}
