//! Physical parameters for simulations

use crate::error::Result;
use crate::physics::values::{DimensionlessNumber, DimensionlessType, Pressure, Velocity};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::geometry::Vector3;
use serde::{Deserialize, Serialize};

/// Physical parameters for simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicalParameters<T: FloatElement + Copy> {
    /// Reynolds number
    pub reynolds_number: DimensionlessNumber<T>,
    /// Reference velocity
    pub reference_velocity: Velocity<T>,
    /// Reference pressure
    pub reference_pressure: Pressure<T>,
    /// Gravity vector
    pub gravity: Vector3<T>,
    /// Time step
    pub time_step: T,
    /// Maximum simulation time
    pub max_time: T,
}

impl<T: FloatElement + Copy + RealField> Default for PhysicalParameters<T> {
    fn default() -> Self {
        Self {
            reynolds_number: DimensionlessNumber::new(
                <T as FloatElement>::from_f64(100.0),
                DimensionlessType::Reynolds,
            )
            .expect("invariant: positive default Reynolds number is valid"),
            reference_velocity: Velocity::from_components(
                <T as NumericElement>::ONE,
                <T as NumericElement>::ZERO,
                <T as NumericElement>::ZERO,
            ),
            reference_pressure: Pressure::from_pascals(<T as FloatElement>::from_f64(101_325.0))
                .expect("invariant: positive default reference pressure is valid"),
            gravity: Vector3::new(
                <T as NumericElement>::ZERO,
                <T as FloatElement>::from_f64(-9.81),
                <T as NumericElement>::ZERO,
            ),
            time_step: <T as FloatElement>::from_f64(0.001),
            max_time: <T as FloatElement>::from_f64(10.0),
        }
    }
}

impl<T: FloatElement + Copy + RealField> PhysicalParameters<T> {
    /// Create parameters with Reynolds number
    ///
    /// # Errors
    ///
    /// Returns an error if `reynolds` violates the dimensionless-number invariant.
    pub fn with_reynolds(reynolds: T) -> Result<Self> {
        Ok(Self {
            reynolds_number: DimensionlessNumber::new(reynolds, DimensionlessType::Reynolds)?,
            ..Default::default()
        })
    }

    /// Set reference velocity
    #[must_use]
    pub fn with_velocity(mut self, velocity: Velocity<T>) -> Self {
        self.reference_velocity = velocity;
        self
    }

    /// Set reference pressure
    #[must_use]
    pub fn with_pressure(mut self, pressure: Pressure<T>) -> Self {
        self.reference_pressure = pressure;
        self
    }

    /// Set gravity vector
    #[must_use]
    pub fn with_gravity(mut self, gravity: Vector3<T>) -> Self {
        self.gravity = gravity;
        self
    }

    /// Set time parameters
    #[must_use]
    pub fn with_time(mut self, dt: T, max_time: T) -> Self {
        self.time_step = dt;
        self.max_time = max_time;
        self
    }

    /// Get CFL number for stability check
    pub fn cfl_number(&self, dx: T) -> T {
        let u_max = self.reference_velocity.magnitude();
        u_max * self.time_step / dx
    }

    /// Check if parameters are stable
    pub fn is_stable(&self, dx: T) -> bool {
        let cfl = self.cfl_number(dx);
        cfl < <T as NumericElement>::ONE
    }
}
