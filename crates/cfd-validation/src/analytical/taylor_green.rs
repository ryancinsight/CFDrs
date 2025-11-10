//! Taylor-Green vortex - decaying vortex solution

use super::AnalyticalSolution;
use nalgebra::{RealField, Vector3};
use cfd_core::conversion::SafeFromF64;
use num_traits::FromPrimitive;
use std::f64::consts::PI;

/// Taylor-Green vortex analytical solution
///
/// Represents a decaying vortex flow that is an exact solution to the
/// incompressible Navier-Stokes equations in 2D/3D periodic domains.
pub struct TaylorGreenVortex<T: RealField + Copy> {
    /// Characteristic length scale
    pub length_scale: T,
    /// Characteristic velocity scale
    pub velocity_scale: T,
    /// Kinematic viscosity
    pub viscosity: T,
    /// Density (for pressure calculation)
    pub density: T,
    /// Whether to use 3D version
    pub is_3d: bool,
}

impl<T: RealField + Copy + FromPrimitive> TaylorGreenVortex<T> {
    /// Create Taylor-Green vortex solution
    pub fn create(
        length_scale: T,
        velocity_scale: T,
        viscosity: T,
        density: T,
        is_3d: bool,
    ) -> Self {
        Self {
            length_scale,
            velocity_scale,
            viscosity,
            density,
            is_3d,
        }
    }

    /// Create 2D Taylor-Green vortex
    pub fn create_2d(length_scale: T, velocity_scale: T, viscosity: T) -> Self {
        Self {
            length_scale,
            velocity_scale,
            viscosity,
            density: T::one(),
            is_3d: false,
        }
    }

    /// Create 3D Taylor-Green vortex
    pub fn create_3d(length_scale: T, velocity_scale: T, viscosity: T) -> Self {
        Self {
            length_scale,
            velocity_scale,
            viscosity,
            density: T::one(),
            is_3d: true,
        }
    }

    /// Get Reynolds number
    pub fn reynolds_number(&self) -> T {
        self.velocity_scale * self.length_scale / self.viscosity
    }

    /// Get the decay rate
    pub fn decay_rate(&self) -> T {
        let pi = <T as SafeFromF64>::try_from_f64(PI)
            .unwrap_or(T::from_f64_or_one(3.14159));
        let factor = if self.is_3d {
            T::from_f64_or_one(3.0)
        } else {
            T::from_f64_or_one(2.0)
        };

        factor * self.viscosity * pi * pi / (self.length_scale * self.length_scale)
    }

    /// Get kinetic energy at time t
    pub fn kinetic_energy(&self, t: T) -> T {
        let initial_energy = if self.is_3d {
            // E₀ = (1/16) * ρ * U² * L³ for 3D
            let factor = <T as SafeFromF64>::try_from_f64(1.0 / 16.0)
                .unwrap_or(T::from_f64_or_one(0.0625));
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * self.length_scale
                * self.length_scale
                * self.length_scale
        } else {
            // E₀ = (1/4) * ρ * U² * L² for 2D
            let factor = T::try_from_f64(0.25)
                .unwrap_or(T::one() / T::from_f64_or_one(4.0));
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * self.length_scale
                * self.length_scale
        };

        let decay = (-(T::from_f64_or_one(2.0)) * self.decay_rate() * t).exp();
        initial_energy * decay
    }

    /// Get enstrophy (vorticity squared) at time t
    pub fn enstrophy(&self, t: T) -> T {
        let pi = <T as SafeFromF64>::try_from_f64(PI)
            .unwrap_or(T::from_f64_or_one(3.14159));
        let initial_enstrophy = self.velocity_scale * self.velocity_scale * pi * pi
            / (self.length_scale * self.length_scale);

        let decay = (-(T::from_f64_or_one(2.0)) * self.decay_rate() * t).exp();
        initial_enstrophy * decay
    }
}

impl<T: RealField + Copy + FromPrimitive> AnalyticalSolution<T> for TaylorGreenVortex<T> {
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        let pi = <T as SafeFromF64>::try_from_f64(PI)
            .unwrap_or(T::from_f64_or_one(3.14159));
        let decay = (-self.decay_rate() * t).exp();

        // Normalize coordinates
        let kx = pi * x / self.length_scale;
        let ky = pi * y / self.length_scale;

        if self.is_3d {
            let kz = pi * z / self.length_scale;

            // 3D Taylor-Green vortex
            let u = self.velocity_scale * kx.sin() * ky.cos() * kz.cos() * decay;
            let v = -self.velocity_scale * kx.cos() * ky.sin() * kz.cos() * decay;
            let w = T::zero(); // For the standard case

            Vector3::new(u, v, w)
        } else {
            // 2D Taylor-Green vortex
            let u = self.velocity_scale * kx.cos() * ky.sin() * decay;
            let v = -self.velocity_scale * kx.sin() * ky.cos() * decay;

            Vector3::new(u, v, T::zero())
        }
    }

    fn pressure(&self, x: T, y: T, z: T, t: T) -> T {
        let pi = <T as SafeFromF64>::try_from_f64(PI)
            .unwrap_or(T::from_f64_or_one(3.14159));
        let decay = (-T::from_f64_or_one(2.0) * self.decay_rate() * t).exp();

        // Normalize coordinates
        let kx = pi * x / self.length_scale;
        let ky = pi * y / self.length_scale;

        if self.is_3d {
            let kz = pi * z / self.length_scale;
            let factor = <T as SafeFromF64>::try_from_f64(1.0 / 16.0)
                .unwrap_or(T::from_f64_or_one(0.0625));

            // p = ρU²/16 * (cos(2kx) + cos(2ky)) * (cos(2kz) + 2) * exp(-2νk²t)
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * (((T::from_f64(2.0).unwrap_or(T::one() + T::one())) * kx).cos()
                    + ((T::from_f64_or_one(2.0)) * ky).cos())
                * (((T::from_f64_or_one(2.0)) * kz).cos()
                    + T::from_f64_or_one(2.0))
                * decay
        } else {
            let factor = T::try_from_f64(0.25)
                .unwrap_or(T::one() / T::from_f64_or_one(4.0));

            // p = -ρU²/4 * (cos(2kx) + cos(2ky)) * exp(-2νk²t)
            -factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * (((T::from_f64_or_one(2.0)) * kx).cos()
                    + ((T::from_f64_or_one(2.0)) * ky).cos())
                * decay
        }
    }

    fn name(&self) -> &str {
        if self.is_3d {
            "3D Taylor-Green Vortex"
        } else {
            "2D Taylor-Green Vortex"
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let two_pi_l = <T as SafeFromF64>::try_from_f64(2.0 * PI)
            .unwrap_or(T::from_f64_or_one(6.28318)) * self.length_scale;
        [
            T::zero(),
            two_pi_l, // x: [0, 2πL]
            T::zero(),
            two_pi_l, // y: [0, 2πL]
            T::zero(),
            if self.is_3d { two_pi_l } else { T::zero() }, // z
        ]
    }

    fn length_scale(&self) -> T {
        self.length_scale
    }

    fn velocity_scale(&self) -> T {
        self.velocity_scale
    }
}
