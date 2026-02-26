//! Lagrangian point representation for IBM
//!
//! # Theorem — Peskin Regularised Delta Properties (Peskin 2002)
//!
//! The discrete delta function $\delta_h$ used for Eulerian–Lagrangian
//! interpolation must satisfy:
//!
//! 1. **Partition of unity:** $\sum_j \delta_h(x - x_j) h = 1$
//! 2. **First moment:** $\sum_j x_j \, \delta_h(x - x_j) h = x$
//! 3. **Compact support:** $\delta_h(r) = 0$ for $|r| \geq r_{\text{cut}}$
//!
//! These properties ensure that interpolation and spreading operations
//! conserve zeroth and first moments (force and torque).
//!
//! **Reference:** Peskin, C.S., "The immersed boundary method",
//! Acta Numerica 11, 2002, pp. 479–517.

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Lagrangian point representing an immersed boundary
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LagrangianPoint<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Position of the Lagrangian point
    pub position: Vector3<T>,
    /// Velocity at the Lagrangian point
    pub velocity: Vector3<T>,
    /// Force at the Lagrangian point
    pub force: Vector3<T>,
    /// Area/volume associated with this point
    pub weight: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> LagrangianPoint<T> {
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
