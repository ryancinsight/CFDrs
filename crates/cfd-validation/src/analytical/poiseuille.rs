//! Poiseuille flow - laminar flow between parallel plates or in a pipe

use super::AnalyticalSolution;
use crate::scalar;
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

/// Poiseuille flow analytical solution
///
/// Represents laminar flow driven by a pressure gradient in a channel or pipe.
/// For parallel plates: u(y) = (`u_max)(1` - (y/h)²)
/// For circular pipe: u(r) = (`u_max)(1` - (r/R)²)
pub struct PoiseuilleFlow<T: RealField + Copy> {
    /// Maximum velocity at centerline
    pub u_max: T,
    /// Channel half-width (plates) or pipe radius
    pub characteristic_length: T,
    /// Pressure gradient magnitude (positive value representing pressure drop per unit length)
    /// For flow in positive x direction, the actual dp/dx is negative
    pub pressure_gradient: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Geometry type
    pub geometry: PoiseuilleGeometry,
}

impl<T: RealField + Copy + FloatElement> PoiseuilleFlow<T> {
    /// Create Poiseuille flow solution
    pub fn create(
        u_max: T,
        characteristic_length: T,
        pressure_gradient: T,
        viscosity: T,
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
        pressure_gradient: T,
        characteristic_length: T,
        viscosity: T,
        geometry: PoiseuilleGeometry,
    ) -> T {
        let factor = match geometry {
            PoiseuilleGeometry::Plates => scalar::from_f64::<T>(2.0),
            PoiseuilleGeometry::Pipe => scalar::from_f64::<T>(4.0),
        };

        scalar::abs(pressure_gradient) * characteristic_length * characteristic_length
            / (factor * viscosity)
    }

    /// Get the flow rate per unit width (plates) or total flow rate (pipe)
    pub fn flow_rate(&self) -> T {
        match self.geometry {
            PoiseuilleGeometry::Plates => {
                // Q = (2/3) * u_max * h
                let two_thirds = scalar::from_f64::<T>(2.0 / 3.0);
                two_thirds * self.u_max * self.characteristic_length
            }
            PoiseuilleGeometry::Pipe => {
                // Q = (π/2) * u_max * R²
                let pi_half = scalar::from_f64::<T>(std::f64::consts::PI / 2.0);
                pi_half * self.u_max * self.characteristic_length * self.characteristic_length
            }
        }
    }

    /// Get Reynolds number
    pub fn reynolds_number(&self, density: T) -> T {
        let characteristic_velocity = self.u_max;
        let two = scalar::from_f64::<T>(2.0);
        let characteristic_length = self.characteristic_length * two; // Diameter for pipe, full width for plates

        density * characteristic_velocity * characteristic_length / self.viscosity
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for PoiseuilleFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        let velocity = match self.geometry {
            PoiseuilleGeometry::Plates => {
                // u(y) = u_max * (1 - (y/h)²)
                let normalized_y = y / self.characteristic_length;
                self.u_max * (scalar::one::<T>() - normalized_y * normalized_y)
            }
            PoiseuilleGeometry::Pipe => {
                // u(r) = u_max * (1 - (r/R)²) where r = sqrt(y² + z²)
                let r = scalar::sqrt(y * y + _z * _z);
                let normalized_r = r / self.characteristic_length;
                if normalized_r <= scalar::one::<T>() {
                    self.u_max * (scalar::one::<T>() - normalized_r * normalized_r)
                } else {
                    scalar::zero::<T>()
                }
            }
        };

        Vector3::new(velocity, scalar::zero::<T>(), scalar::zero::<T>())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop: p(x) = p0 - (dp/dx) * x
        -self.pressure_gradient * x
    }

    fn name(&self) -> &str {
        match self.geometry {
            PoiseuilleGeometry::Plates => "Poiseuille Flow (Parallel Plates)",
            PoiseuilleGeometry::Pipe => "Poiseuille Flow (Pipe)",
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let large = scalar::from_f64::<T>(1000.0);
        match self.geometry {
            PoiseuilleGeometry::Plates => [
                scalar::zero::<T>(),
                large, // x: [0, L]
                -self.characteristic_length,
                self.characteristic_length, // y: [-h, h]
                -large,
                large, // z: arbitrary
            ],
            PoiseuilleGeometry::Pipe => [
                scalar::zero::<T>(),
                large, // x: [0, L]
                -self.characteristic_length,
                self.characteristic_length, // y: [-R, R]
                -self.characteristic_length,
                self.characteristic_length, // z: [-R, R]
            ],
        }
    }

    fn length_scale(&self) -> T {
        self.characteristic_length
    }

    fn velocity_scale(&self) -> T {
        self.u_max
    }
}
