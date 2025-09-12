//! Velocity value object

use crate::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::fmt;

// Velocity conversion constants
const MS_TO_KMH: f64 = 3.6;
const MS_TO_MPH: f64 = 2.23694;
const MS_TO_KNOT: f64 = 1.94384;
const MS_TO_FTS: f64 = 3.28084;

// Mach number thresholds
const SUBSONIC_MACH_LIMIT: f64 = 0.8;
const SUPERSONIC_MACH_LIMIT: f64 = 1.2;

/// Velocity vector with magnitude and direction
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Velocity<T: RealField + Copy> {
    /// Velocity components in m/s (SI units)
    components: Vector3<T>,
}

impl<T: RealField + Copy + FromPrimitive> Velocity<T> {
    /// Create zero velocity
    #[must_use]
    pub fn zero() -> Self {
        Self {
            components: Vector3::zeros(),
        }
    }

    /// Create velocity from components in m/s
    pub fn from_components(x: T, y: T, z: T) -> Self {
        Self {
            components: Vector3::new(x, y, z),
        }
    }

    /// Create velocity from vector
    pub fn from_vector(vec: Vector3<T>) -> Self {
        Self { components: vec }
    }

    /// Create uniform velocity in x-direction
    pub fn uniform_x(speed: T) -> Self {
        Self::from_components(speed, T::zero(), T::zero())
    }

    /// Get velocity components
    pub fn components(&self) -> &Vector3<T> {
        &self.components
    }

    /// Get velocity magnitude
    pub fn magnitude(&self) -> T {
        self.components.norm()
    }

    /// Get velocity direction (unit vector)
    pub fn direction(&self) -> Option<Vector3<T>> {
        let mag = self.magnitude();
        if mag > T::zero() {
            Some(self.components / mag)
        } else {
            None
        }
    }

    /// Get x-component
    pub fn x(&self) -> T {
        self.components.x
    }

    /// Get y-component
    pub fn y(&self) -> T {
        self.components.y
    }

    /// Get z-component
    pub fn z(&self) -> T {
        self.components.z
    }

    /// Get magnitude in km/h
    pub fn kmh(&self) -> T {
        self.magnitude() * T::from_f64(MS_TO_KMH).unwrap_or_else(T::zero)
    }

    /// Get magnitude in mph
    pub fn mph(&self) -> T {
        self.magnitude() * T::from_f64(MS_TO_MPH).unwrap_or_else(T::zero)
    }

    /// Get magnitude in knots
    pub fn knots(&self) -> T {
        self.magnitude() * T::from_f64(MS_TO_KNOT).unwrap_or_else(T::zero)
    }

    /// Get magnitude in ft/s
    pub fn fps(&self) -> T {
        self.magnitude() * T::from_f64(MS_TO_FTS).unwrap_or_else(T::zero)
    }

    /// Check if velocity is subsonic (Mach < `SUBSONIC_MACH_LIMIT`)
    pub fn is_subsonic(&self, speed_of_sound: T) -> bool {
        let mach_limit = T::from_f64(SUBSONIC_MACH_LIMIT).unwrap_or_else(T::zero);
        self.magnitude() < mach_limit * speed_of_sound
    }

    /// Check if velocity is supersonic (Mach > `SUPERSONIC_MACH_LIMIT`)
    pub fn is_supersonic(&self, speed_of_sound: T) -> bool {
        let mach_limit = T::from_f64(SUPERSONIC_MACH_LIMIT).unwrap_or_else(T::zero);
        self.magnitude() > mach_limit * speed_of_sound
    }

    /// Calculate Mach number
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn mach_number(&self, speed_of_sound: T) -> Result<T> {
        if speed_of_sound <= T::zero() {
            return Err(crate::error::Error::InvalidConfiguration(
                "Speed of sound must be positive".into(),
            ));
        }
        Ok(self.magnitude() / speed_of_sound)
    }
}

impl<T: RealField + Copy + fmt::Display> fmt::Display for Velocity<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({}, {}, {}) m/s",
            self.components.x, self.components.y, self.components.z
        )
    }
}
