//! Wall boundary condition types
//!
//! Reference: Versteeg & Malalasekera (2007), Ch. 9

use nalgebra::{RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Wall boundary condition types
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum WallType<T: RealField + Copy> {
    /// No-slip wall (u = 0)
    NoSlip,

    /// Slip wall (u·n = 0)
    Slip,

    /// Moving wall with prescribed velocity
    Moving {
        /// Wall velocity vector
        velocity: Vector3<T>,
    },

    /// Rotating wall
    Rotating {
        /// Angular velocity vector [rad/s]
        omega: Vector3<T>,
        /// Center of rotation
        center: Vector3<T>,
    },
}

impl<T: RealField + Copy> WallType<T> {
    /// Create no-slip wall
    #[must_use]
    pub const fn no_slip() -> Self {
        Self::NoSlip
    }

    /// Create slip wall
    #[must_use]
    pub const fn slip() -> Self {
        Self::Slip
    }

    /// Create moving wall
    pub fn moving(velocity: Vector3<T>) -> Self {
        Self::Moving { velocity }
    }

    /// Create rotating wall
    pub fn rotating(omega: Vector3<T>, center: Vector3<T>) -> Self {
        Self::Rotating { omega, center }
    }
}
