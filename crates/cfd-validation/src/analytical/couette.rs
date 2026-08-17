//! Couette flow - shear-driven flow between parallel plates

use super::AnalyticalSolution;
use crate::scalar;
use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Length, MassDensity, Pressure, PressureGradient,
    ReciprocalTime, Velocity,
};
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;

/// Couette flow analytical solution
///
/// Represents laminar flow between parallel plates where one plate moves
/// with velocity U relative to the other. Can include pressure gradient effects.
pub struct CouetteFlow<T: RealField + Copy> {
    /// Velocity of the moving wall
    pub wall_velocity: Velocity<T>,
    /// Gap height between plates
    pub gap_height: Length<T>,
    /// Pressure gradient (optional)
    pub pressure_gradient: PressureGradient<T>,
    /// Dynamic viscosity
    pub viscosity: DynamicViscosity<T>,
}

impl<T: RealField + Copy + FloatElement> CouetteFlow<T> {
    /// Create Couette flow solution
    pub fn create(
        wall_velocity: Velocity<T>,
        gap_height: Length<T>,
        pressure_gradient: PressureGradient<T>,
        viscosity: DynamicViscosity<T>,
    ) -> Self {
        Self {
            wall_velocity,
            gap_height,
            pressure_gradient,
            viscosity,
        }
    }

    /// Create pure Couette flow (no pressure gradient)
    pub fn pure(wall_velocity: Velocity<T>, gap_height: Length<T>) -> Self {
        Self {
            wall_velocity,
            gap_height,
            pressure_gradient: PressureGradient::from_base(scalar::zero::<T>()),
            viscosity: DynamicViscosity::from_base(scalar::one::<T>()),
        }
    }

    /// Get the shear rate
    pub fn shear_rate(&self) -> ReciprocalTime<T> {
        ReciprocalTime::from_base(self.wall_velocity.into_base() / self.gap_height.into_base())
    }

    /// Get the wall shear stress
    pub fn wall_shear_stress(&self) -> Pressure<T> {
        let base_shear = self.viscosity.into_base() * self.wall_velocity.into_base()
            / self.gap_height.into_base();
        let pressure_contribution = self.pressure_gradient.into_base()
            * self.gap_height.into_base()
            / <T as FloatElement>::from_f64(2.0);
        Pressure::from_base(base_shear + pressure_contribution)
    }

    /// Get Reynolds number based on gap height
    pub fn reynolds_number(&self, density: MassDensity<T>) -> Dimensionless<T> {
        Dimensionless::from_base(
            density.into_base() * self.wall_velocity.into_base() * self.gap_height.into_base()
                / self.viscosity.into_base(),
        )
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for CouetteFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        // Normalize y coordinate: η = y/h where y ∈ [0, h]
        let gap_height = self.gap_height.into_base();
        let wall_velocity = self.wall_velocity.into_base();
        let pressure_gradient = self.pressure_gradient.into_base();
        let viscosity = self.viscosity.into_base();
        let eta = y / gap_height;

        // Couette-Poiseuille flow: u(y) = U*y/h + (1/2μ)(dp/dx)(y)(y-h)
        let couette_part = wall_velocity * eta;

        let poiseuille_part = if pressure_gradient == scalar::zero::<T>() {
            scalar::zero::<T>()
        } else {
            // Plane Poiseuille contribution (Versteeg & Malalasekera):
            // u_p(y) = -(1/(2μ)) (dp/dx) y (h - y)
            let two = <T as FloatElement>::from_f64(2.0);
            let factor = -pressure_gradient / (two * viscosity);
            factor * y * (gap_height - y)
        };

        let u = couette_part + poiseuille_part;

        Vector3::new(u, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop if pressure gradient exists
        -self.pressure_gradient.into_base() * x
    }

    fn name(&self) -> &str {
        if self.pressure_gradient.into_base() == scalar::zero::<T>() {
            "Pure Couette Flow"
        } else {
            "Couette-Poiseuille Flow"
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let large = <T as FloatElement>::from_f64(1000.0);
        [
            scalar::zero::<T>(),
            large, // x: [0, L]
            scalar::zero::<T>(),
            self.gap_height.into_base(), // y: [0, h]
            -large,
            large, // z: arbitrary
        ]
    }

    fn length_scale(&self) -> T {
        self.gap_height.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.wall_velocity.into_base()
    }
}
