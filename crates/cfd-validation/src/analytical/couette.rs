//! Couette flow - shear-driven flow between parallel plates

use super::AnalyticalSolution;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Couette flow analytical solution
///
/// Represents laminar flow between parallel plates where one plate moves
/// with velocity U relative to the other. Can include pressure gradient effects.
pub struct CouetteFlow<T: RealField + Copy> {
    /// Velocity of the moving wall
    pub wall_velocity: T,
    /// Gap height between plates
    pub gap_height: T,
    /// Pressure gradient (optional)
    pub pressure_gradient: T,
    /// Dynamic viscosity
    pub viscosity: T,
}

impl<T: RealField + Copy + FromPrimitive> CouetteFlow<T> {
    /// Create Couette flow solution
    pub fn create(wall_velocity: T, gap_height: T, pressure_gradient: T, viscosity: T) -> Self {
        Self {
            wall_velocity,
            gap_height,
            pressure_gradient,
            viscosity,
        }
    }

    /// Create pure Couette flow (no pressure gradient)
    pub fn pure(wall_velocity: T, gap_height: T) -> Self {
        Self {
            wall_velocity,
            gap_height,
            pressure_gradient: T::zero(),
            viscosity: T::one(), // Normalized
        }
    }

    /// Get the shear rate
    pub fn shear_rate(&self) -> T {
        self.wall_velocity / self.gap_height
    }

    /// Get the wall shear stress
    pub fn wall_shear_stress(&self) -> T {
        let base_shear = self.viscosity * self.wall_velocity / self.gap_height;
        let pressure_contribution = self.pressure_gradient * self.gap_height
            / T::from_f64(2.0).unwrap_or(T::one() + T::one());
        base_shear + pressure_contribution
    }

    /// Get Reynolds number based on gap height
    pub fn reynolds_number(&self, density: T) -> T {
        density * self.wall_velocity * self.gap_height / self.viscosity
    }
}

impl<T: RealField + Copy + FromPrimitive> AnalyticalSolution<T> for CouetteFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        // Normalize y coordinate: η = y/h where y ∈ [0, h]
        let eta = y / self.gap_height;

        // Couette-Poiseuille flow: u(y) = U*y/h + (1/2μ)(dp/dx)(y)(y-h)
        let couette_part = self.wall_velocity * eta;

        let poiseuille_part = if self.pressure_gradient == T::zero() {
            T::zero()
        } else {
            // Plane Poiseuille contribution (Versteeg & Malalasekera):
            // u_p(y) = -(1/(2μ)) (dp/dx) y (h - y)
            let two = T::from_f64(2.0).unwrap_or(T::one() + T::one());
            let factor = -self.pressure_gradient / (two * self.viscosity);
            factor * y * (self.gap_height - y)
        };

        let u = couette_part + poiseuille_part;

        Vector3::new(u, T::zero(), T::zero())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop if pressure gradient exists
        -self.pressure_gradient * x
    }

    fn name(&self) -> &str {
        if self.pressure_gradient == T::zero() {
            "Pure Couette Flow"
        } else {
            "Couette-Poiseuille Flow"
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let large = T::from_f64(1000.0).unwrap_or(T::from_f64(100.0).unwrap_or(T::one()));
        [
            T::zero(),
            large, // x: [0, L]
            T::zero(),
            self.gap_height, // y: [0, h]
            -large,
            large, // z: arbitrary
        ]
    }

    fn length_scale(&self) -> T {
        self.gap_height
    }

    fn velocity_scale(&self) -> T {
        self.wall_velocity
    }
}
