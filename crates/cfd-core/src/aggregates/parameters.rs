//! Physical parameters for simulations

use crate::values::{DimensionlessNumber, DimensionlessType, Pressure, Velocity};
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// Physical parameters for simulation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhysicalParameters<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive> Default for PhysicalParameters<T> {
    fn default() -> Self {
        Self {
            reynolds_number: DimensionlessNumber::new(
                T::from_f64(100.0).unwrap_or_else(T::one),
                DimensionlessType::Reynolds,
            )
            .unwrap_or_else(|_| {
                DimensionlessNumber::new(T::one(), DimensionlessType::Reynolds).unwrap()
            }),
            reference_velocity: Velocity::from_components(T::one(), T::zero(), T::zero()),
            reference_pressure: Pressure::from_pascals(
                T::from_f64(101_325.0).unwrap_or_else(T::one),
            )
            .unwrap_or_else(|_| Pressure::zero()),
            gravity: Vector3::new(
                T::zero(),
                T::from_f64(-9.81).unwrap_or_else(T::zero),
                T::zero(),
            ),
            time_step: T::from_f64(0.001)
                .unwrap_or_else(|| T::from_f64(0.001).unwrap_or_else(T::one)),
            max_time: T::from_f64(10.0).unwrap_or_else(|| T::from_f64(10.0).unwrap_or_else(T::one)),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> PhysicalParameters<T> {
    /// Create parameters with Reynolds number
    ///
    /// # Panics
    /// May panic if the default Reynolds number creation fails due to invalid numeric type conversion
    pub fn with_reynolds(reynolds: T) -> Self {
        Self {
            reynolds_number: DimensionlessNumber::new(reynolds, DimensionlessType::Reynolds)
                .unwrap_or_else(|_| {
                    DimensionlessNumber::new(T::one(), DimensionlessType::Reynolds).unwrap()
                }),
            ..Default::default()
        }
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
        cfl < T::one()
    }
}
