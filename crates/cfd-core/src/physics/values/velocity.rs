//! Velocity value object with an aequitas-backed unit layer.
//!
//! The direction-bearing core stays a `leto` 3-vector in m/s (geometry, not
//! units); every scalar magnitude conversion (km/h, mph, knots, ft/s) goes
//! through the aequitas velocity catalog rather than local constants.

use super::dimensionless::DimensionlessNumber;
use crate::error::Result;
use crate::physics::constants::physics::dimensionless::mach::{TRANSONIC_LOWER, TRANSONIC_UPPER};
use aequitas::systems::si::quantities::Velocity as VelocityQuantity;
use aequitas::systems::si::units::{
    FootPerSecond, KilometerPerHour, Knot, MeterPerSecond, MilePerHour,
};
use eunomia::{FloatElement, NumericElement, RealField, UnitScalar};
use leto::geometry::Vector3;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Velocity vector with magnitude and direction.
///
/// Components are stored in m/s. Construct from raw components at the solver
/// boundary ([`from_components`](Self::from_components)); read magnitudes in
/// any catalog unit through the `kmh`/`mph`/`knots`/`fps` accessors.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Velocity<T: FloatElement + Copy> {
    /// Velocity components in m/s (SI units)
    components: Vector3<T>,
}

impl<T: FloatElement + Copy + RealField> Velocity<T> {
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
        Self::from_components(
            speed,
            <T as NumericElement>::ZERO,
            <T as NumericElement>::ZERO,
        )
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
        if mag > <T as NumericElement>::ZERO {
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

    /// Check if velocity is subsonic (Mach below the central transonic lower
    /// edge, 0.8).
    pub fn is_subsonic(&self, speed_of_sound: T) -> bool {
        let mach_limit = <T as FloatElement>::from_f64(TRANSONIC_LOWER);
        self.magnitude() < mach_limit * speed_of_sound
    }

    /// Check if velocity is supersonic (Mach above the central transonic upper
    /// edge, 1.2).
    pub fn is_supersonic(&self, speed_of_sound: T) -> bool {
        let mach_limit = <T as FloatElement>::from_f64(TRANSONIC_UPPER);
        self.magnitude() > mach_limit * speed_of_sound
    }

    /// Calculate Mach number.
    ///
    /// The ratio uses the single [`DimensionlessNumber::mach`] formula; this
    /// method only supplies the vector magnitude. A non-negative magnitude
    /// over a positive speed of sound can never trip the taxonomy's
    /// non-negativity guard, so the error paths coincide exactly.
    ///
    /// # Errors
    ///
    /// Returns an error if the input value is not finite or is outside valid range.
    pub fn mach_number(&self, speed_of_sound: T) -> Result<T> {
        DimensionlessNumber::mach(self.magnitude(), speed_of_sound).map(|mach| mach.value())
    }
}

impl<T: FloatElement + Copy + RealField + UnitScalar> Velocity<T> {
    /// Magnitude as an aequitas velocity quantity (m/s coherent).
    #[must_use]
    pub fn quantity(&self) -> VelocityQuantity<T> {
        VelocityQuantity::from_base(self.magnitude())
    }

    /// Get magnitude in km/h.
    pub fn kmh(&self) -> T {
        self.quantity().in_unit::<KilometerPerHour>()
    }

    /// Get magnitude in mph.
    pub fn mph(&self) -> T {
        self.quantity().in_unit::<MilePerHour>()
    }

    /// Get magnitude in knots.
    pub fn knots(&self) -> T {
        self.quantity().in_unit::<Knot>()
    }

    /// Get magnitude in ft/s.
    pub fn fps(&self) -> T {
        self.quantity().in_unit::<FootPerSecond>()
    }

    /// Get magnitude in m/s.
    pub fn meters_per_second(&self) -> T {
        self.quantity().in_unit::<MeterPerSecond>()
    }
}

impl<T: FloatElement + Copy + fmt::Display> fmt::Display for Velocity<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "({}, {}, {}) m/s",
            self.components.x, self.components.y, self.components.z
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmh_mph_knots_fps_match_definitions() {
        // 10 m/s parade case: 36 km/h; 22.36936... mph; 19.43844... kn.
        let v = Velocity::<f64>::uniform_x(10.0);
        assert!((v.kmh() - 36.0).abs() < 1e-9);
        assert!((v.mph() - 22.369_362_920_8).abs() < 1e-9);
        assert!((v.knots() - 19.438_444_924_4).abs() < 1e-9);
        assert!((v.fps() - 32.808_398_950_1).abs() < 1e-9);
        assert!((v.meters_per_second() - 10.0).abs() < 1e-12);
    }

    #[test]
    fn mach_regimes_classify() {
        let v = Velocity::<f64>::uniform_x(100.0);
        assert!(v.is_subsonic(340.0));
        assert!(!v.is_supersonic(340.0));
        assert!((v.mach_number(340.0).expect("valid") - 100.0 / 340.0).abs() < 1e-12);
        assert!(v.mach_number(0.0).is_err());
    }

    #[test]
    fn mach_number_matches_taxonomy_formula() {
        // The vector method must agree bit-for-bit with the scalar taxonomy
        // constructor it delegates to, including the error path.
        let v = Velocity::<f64>::uniform_x(100.0);
        let taxonomy = DimensionlessNumber::<f64>::mach(100.0, 340.0).expect("valid");
        assert_eq!(
            v.mach_number(340.0).expect("valid").to_bits(),
            taxonomy.value().to_bits()
        );
        assert!(DimensionlessNumber::mach(100.0, 0.0).is_err());
    }
}
