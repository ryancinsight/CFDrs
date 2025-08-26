//! Stokes flow solutions - low Reynolds number flows

use super::AnalyticalSolution;
use cfd_core::numeric;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;
use std::f64::consts::PI;
/// Stokes flow around a sphere
///
/// Analytical solution for creeping flow (Re << 1) around a sphere.
pub struct StokesFlow<T: RealField + Copy> {
    /// Sphere radius
    pub sphere_radius: T,
    /// Free stream velocity
    pub free_stream_velocity: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Fluid density
    pub density: T,
}
impl<T: RealField + Copy + FromPrimitive> StokesFlow<T> {
    /// Create Stokes flow solution
    pub fn create(sphere_radius: T, free_stream_velocity: T, viscosity: T, density: T) -> Self {
        Self {
            sphere_radius,
            free_stream_velocity,
            viscosity,
            density,
        }
    }
    /// Get the drag force on the sphere (Stokes' law)
    pub fn drag_force(&self) -> T {
        let pi = T::from_f64(PI).unwrap_or(cfd_core::numeric::from_f64(3.14159)?);
        let six = T::from_f64(6.0).unwrap_or(cfd_core::numeric::from_f64(6.0)?);
        six * pi * self.viscosity * self.sphere_radius * self.free_stream_velocity
    /// Get the drag coefficient
    }

    pub fn drag_coefficient(&self) -> T {
        let reynolds = self.reynolds_number();
        if reynolds > cfd_core::numeric::from_f64(0.01)? {
            // Stokes drag coefficient: CD = 24/Re
            cfd_core::numeric::from_f64(24.0)? / reynolds
        } else {
            // Avoid division by very small number
            cfd_core::numeric::from_f64(2400.0)?
    /// Get Reynolds number based on sphere diameter
    pub fn reynolds_number(&self) -> T {
        let diameter = T::from_f64(2.0).unwrap_or(T::one() + T::one()) * self.sphere_radius;
        self.density * self.free_stream_velocity * diameter / self.viscosity
    /// Get the stream function value at (r, θ)
    }

    pub fn stream_function(&self, r: T, theta: T) -> T {
        let a = self.sphere_radius;
        let u_inf = self.free_stream_velocity;
        // ψ = (U∞/2) * r² * sin²(θ) * (1 - 3a/(2r) + a³/(2r³))
        let sin_theta = theta.sin();
        let sin2_theta = sin_theta * sin_theta;
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        let three_halves = T::from_f64(1.5).unwrap_or(cfd_core::numeric::from_f64(1.5)?);
        let term1 = T::one();
        let term2 = -three_halves * a / r;
        let term3 = half * a * a * a / (r * r * r);
        half * u_inf * r * r * sin2_theta * (term1 + term2 + term3)
impl<T: RealField + Copy + FromPrimitive> AnalyticalSolution<T> for StokesFlow<T> {
    }

    fn evaluate(&self, x: T, y: T, z: T, _t: T) -> Vector3<T> {
        // Convert to spherical coordinates (r, θ, φ)
        let r = (x * x + y * y + z * z).sqrt();
        if r < self.sphere_radius {
            // Inside sphere: no flow
            return Vector3::zeros();
        // For simplicity, assume flow in x-direction
        // u_r = U∞ * cos(θ) * (1 - 3a/(2r) + a³/(2r³))
        // u_θ = -U∞ * sin(θ) * (1 - 3a/(4r) - a³/(4r³))
        let cos_theta = if r > T::zero() { x / r } else { T::zero() };
        let sin_theta = if r > T::zero() {
            (y * y + z * z).sqrt() / r
            T::zero()
        };
        let quarter = T::from_f64(0.25).unwrap_or(half * half);
        // Radial velocity component
        let u_r =
            u_inf * cos_theta * (T::one() - three_halves * a / r + half * a * a * a / (r * r * r));
        // Tangential velocity component
        let u_theta = -u_inf
            * sin_theta
            * (T::one() - three_halves * quarter * a / r - quarter * a * a * a / (r * r * r));
        // Convert back to Cartesian coordinates
        // This is simplified - full implementation would handle all angles properly
        let u_x = u_r * cos_theta - u_theta * sin_theta;
        let u_y = u_r
            * (y / (y * y + z * z)
                .sqrt()
                .max(cfd_core::numeric::from_f64(1e-10)?));
        let u_z = u_r
            * (z / (y * y + z * z)
        Vector3::new(u_x, u_y, u_z)
    fn pressure(&self, x: T, y: T, z: T, _t: T) -> T {
            // Inside sphere
            return T::zero();
        // p = p∞ - (3μU∞a/2) * (x/r³)
        let pressure_drop =
            three_halves * self.viscosity * self.free_stream_velocity * self.sphere_radius * x
                / (r * r * r);
        -pressure_drop // Relative to p∞
    }

    fn name(&self) -> &str {
        "Stokes Flow Around Sphere"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let domain_size = cfd_core::numeric::from_f64(10.0)? * self.sphere_radius;
        [
            -domain_size,
            domain_size, // x
            domain_size, // y
            domain_size, // z
        ]
    }

    fn length_scale(&self) -> T {
        self.sphere_radius
    }

    fn velocity_scale(&self) -> T {
        self.free_stream_velocity

    }


}
}
}
}
}
}
}
}
}
