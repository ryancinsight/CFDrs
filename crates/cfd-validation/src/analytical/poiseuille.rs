//! Poiseuille flow - laminar flow between parallel plates or in a pipe

use super::AnalyticalSolution;
use crate::scalar;
use aequitas::systems::si::quantities::{
    AreaPerTime, Dimensionless, DynamicViscosity, Length, MassDensity, PressureGradient, Velocity,
    VolumetricFlowRate,
};
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;

/// Poiseuille flow configuration
#[derive(Debug, Clone, Copy)]
pub enum PoiseuilleGeometry {
    /// Flow between parallel plates
    Plates,
    /// Flow in a circular pipe
    Pipe,
}

/// Flow-rate result whose unit follows the selected geometry.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PoiseuilleFlowRate<T> {
    /// Planar flow rate per unit width, in square metres per second.
    PerWidth(AreaPerTime<T>),
    /// Pipe flow rate, in cubic metres per second.
    Volumetric(VolumetricFlowRate<T>),
}

/// Poiseuille flow analytical solution
///
/// Represents laminar flow driven by a pressure gradient in a channel or pipe.
/// For parallel plates: u(y) = (`u_max)(1` - (y/h)²)
/// For circular pipe: u(r) = (`u_max)(1` - (r/R)²)
pub struct PoiseuilleFlow<T: RealField + Copy> {
    /// Maximum velocity at centerline
    pub u_max: Velocity<T>,
    /// Channel half-width (plates) or pipe radius
    pub characteristic_length: Length<T>,
    /// Pressure gradient magnitude (positive value representing pressure drop per unit length)
    /// For flow in positive x direction, the actual dp/dx is negative
    pub pressure_gradient: PressureGradient<T>,
    /// Dynamic viscosity
    pub viscosity: DynamicViscosity<T>,
    /// Geometry type
    pub geometry: PoiseuilleGeometry,
}

impl<T: RealField + Copy + FloatElement> PoiseuilleFlow<T> {
    /// Create Poiseuille flow solution
    pub fn create(
        u_max: Velocity<T>,
        characteristic_length: Length<T>,
        pressure_gradient: PressureGradient<T>,
        viscosity: DynamicViscosity<T>,
        geometry: PoiseuilleGeometry,
    ) -> Self {
        Self {
            u_max,
            characteristic_length,
            pressure_gradient,
            viscosity,
            geometry,
        }
    }

    /// Calculate velocity from pressure gradient
    pub fn velocity_from_pressure_gradient(
        pressure_gradient: PressureGradient<T>,
        characteristic_length: Length<T>,
        viscosity: DynamicViscosity<T>,
        geometry: PoiseuilleGeometry,
    ) -> Velocity<T> {
        let factor = match geometry {
            PoiseuilleGeometry::Plates => <T as FloatElement>::from_f64(2.0),
            PoiseuilleGeometry::Pipe => <T as FloatElement>::from_f64(4.0),
        };

        Velocity::from_base(
            scalar::abs(pressure_gradient.into_base())
                * characteristic_length.into_base()
                * characteristic_length.into_base()
                / (factor * viscosity.into_base()),
        )
    }

    /// Get the geometry-specific flow rate.
    pub fn flow_rate(&self) -> PoiseuilleFlowRate<T> {
        let u_max = self.u_max.into_base();
        let characteristic_length = self.characteristic_length.into_base();
        match self.geometry {
            PoiseuilleGeometry::Plates => {
                // Q/W = ∫_{-h}^{h} u_max (1 − (y/h)²) dy = (4/3)·u_max·h,
                // equivalently ū·2h with ū = (2/3)·u_max over the full gap 2h
                // (characteristic_length is the half-gap h; see u_max above).
                let four_thirds = <T as FloatElement>::from_f64(4.0 / 3.0);
                PoiseuilleFlowRate::PerWidth(AreaPerTime::from_base(
                    four_thirds * u_max * characteristic_length,
                ))
            }
            PoiseuilleGeometry::Pipe => {
                // Q = (π/2) * u_max * R²
                let pi_half = <T as FloatElement>::from_f64(std::f64::consts::PI / 2.0);
                PoiseuilleFlowRate::Volumetric(VolumetricFlowRate::from_base(
                    pi_half * u_max * characteristic_length * characteristic_length,
                ))
            }
        }
    }

    /// Get Reynolds number
    pub fn reynolds_number(&self, density: MassDensity<T>) -> Dimensionless<T> {
        let characteristic_velocity = self.u_max.into_base();
        let two = <T as FloatElement>::from_f64(2.0);
        let characteristic_length = self.characteristic_length.into_base() * two;

        Dimensionless::from_base(
            density.into_base() * characteristic_velocity * characteristic_length
                / self.viscosity.into_base(),
        )
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for PoiseuilleFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let velocity = match self.geometry {
            PoiseuilleGeometry::Plates => {
                // u(y) = u_max * (1 - (y/h)²)
                let normalized_y = y / self.characteristic_length.into_base();
                self.u_max.into_base() * (scalar::one::<T>() - normalized_y * normalized_y)
            }
            PoiseuilleGeometry::Pipe => {
                // u(r) = u_max * (1 - (r/R)²) where r = sqrt(y² + z²)
                let r = scalar::sqrt(y * y + _z * _z);
                let normalized_r = r / self.characteristic_length.into_base();
                if normalized_r <= scalar::one::<T>() {
                    self.u_max.into_base() * (scalar::one::<T>() - normalized_r * normalized_r)
                } else {
                    scalar::zero::<T>()
                }
            }
        };

        Vector3::new(velocity, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop: p(x) = p0 - (dp/dx) * x
        -self.pressure_gradient.into_base() * x
    }

    fn name(&self) -> &str {
        match self.geometry {
            PoiseuilleGeometry::Plates => "Poiseuille Flow (Parallel Plates)",
            PoiseuilleGeometry::Pipe => "Poiseuille Flow (Pipe)",
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let large = <T as FloatElement>::from_f64(1000.0);
        match self.geometry {
            PoiseuilleGeometry::Plates => [
                scalar::zero::<T>(),
                large, // x: [0, L]
                -self.characteristic_length.into_base(),
                self.characteristic_length.into_base(), // y: [-h, h]
                -large,
                large, // z: arbitrary
            ],
            PoiseuilleGeometry::Pipe => [
                scalar::zero::<T>(),
                large, // x: [0, L]
                -self.characteristic_length.into_base(),
                self.characteristic_length.into_base(), // y: [-R, R]
                -self.characteristic_length.into_base(),
                self.characteristic_length.into_base(), // z: [-R, R]
            ],
        }
    }

    fn length_scale(&self) -> T {
        self.characteristic_length.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.u_max.into_base()
    }
}
