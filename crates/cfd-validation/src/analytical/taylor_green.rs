//! Taylor-Green vortex - decaying vortex solution

use super::AnalyticalSolution;
use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;
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

impl<T: RealField + Copy + FloatElement> TaylorGreenVortex<T> {
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
            density: scalar::one::<T>(),
            is_3d: false,
        }
    }

    /// Create 3D Taylor-Green vortex
    pub fn create_3d(length_scale: T, velocity_scale: T, viscosity: T) -> Self {
        Self {
            length_scale,
            velocity_scale,
            viscosity,
            density: scalar::one::<T>(),
            is_3d: true,
        }
    }

    /// Get Reynolds number
    pub fn reynolds_number(&self) -> T {
        self.velocity_scale * self.length_scale / self.viscosity
    }

    /// Get the decay rate
    pub fn decay_rate(&self) -> T {
        let pi = scalar::from_f64::<T>(PI);
        let factor = if self.is_3d {
            scalar::from_f64::<T>(3.0)
        } else {
            scalar::from_f64::<T>(2.0)
        };

        factor * self.viscosity * pi * pi / (self.length_scale * self.length_scale)
    }

    /// Get kinetic energy at time t
    pub fn kinetic_energy(&self, t: T) -> T {
        let initial_energy = if self.is_3d {
            // E₀ = (1/16) * ρ * U² * L³ for 3D
            let factor = scalar::from_f64::<T>(1.0 / 16.0);
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * self.length_scale
                * self.length_scale
                * self.length_scale
        } else {
            // E₀ = (1/4) * ρ * U² * L² for 2D
            let factor = scalar::from_f64::<T>(0.25);
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * self.length_scale
                * self.length_scale
        };

        let decay = scalar::exp(-(scalar::from_f64::<T>(2.0)) * self.decay_rate() * t);
        initial_energy * decay
    }

    /// Get enstrophy (vorticity squared) at time t
    pub fn enstrophy(&self, t: T) -> T {
        let pi = scalar::from_f64::<T>(PI);
        let initial_enstrophy = self.velocity_scale * self.velocity_scale * pi * pi
            / (self.length_scale * self.length_scale);

        let decay = scalar::exp(-(scalar::from_f64::<T>(2.0)) * self.decay_rate() * t);
        initial_enstrophy * decay
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for TaylorGreenVortex<T> {
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        let pi = scalar::from_f64::<T>(PI);
        let decay = scalar::exp(-self.decay_rate() * t);

        // Normalize coordinates
        let kx = pi * x / self.length_scale;
        let ky = pi * y / self.length_scale;

        if self.is_3d {
            let kz = pi * z / self.length_scale;

            // 3D Taylor-Green vortex
            let u =
                self.velocity_scale * scalar::sin(kx) * scalar::cos(ky) * scalar::cos(kz) * decay;
            let v =
                -self.velocity_scale * scalar::cos(kx) * scalar::sin(ky) * scalar::cos(kz) * decay;
            let w = scalar::zero::<T>(); // For the standard case

            Vector3::new(u, v, w)
        } else {
            // 2D Taylor-Green vortex
            let u = self.velocity_scale * scalar::cos(kx) * scalar::sin(ky) * decay;
            let v = -self.velocity_scale * scalar::sin(kx) * scalar::cos(ky) * decay;

            Vector3::new(u, v, scalar::zero::<T>())
        }
    }

    fn pressure(&self, x: T, y: T, z: T, t: T) -> T {
        let pi = scalar::from_f64::<T>(PI);
        let decay = scalar::exp(-scalar::from_f64::<T>(2.0) * self.decay_rate() * t);

        // Normalize coordinates
        let kx = pi * x / self.length_scale;
        let ky = pi * y / self.length_scale;

        if self.is_3d {
            let kz = pi * z / self.length_scale;
            let factor = scalar::from_f64::<T>(1.0 / 16.0);
            let two = scalar::from_f64::<T>(2.0);

            // p = ρU²/16 * (cos(2kx) + cos(2ky)) * (cos(2kz) + 2) * exp(-2νk²t)
            factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * (scalar::cos(two * kx) + scalar::cos(two * ky))
                * (scalar::cos(two * kz) + two)
                * decay
        } else {
            let factor = scalar::from_f64::<T>(0.25);
            let two = scalar::from_f64::<T>(2.0);

            // p = -ρU²/4 * (cos(2kx) + cos(2ky)) * exp(-2νk²t)
            -factor
                * self.density
                * self.velocity_scale
                * self.velocity_scale
                * (scalar::cos(two * kx) + scalar::cos(two * ky))
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
        let two_pi_l = scalar::from_f64::<T>(2.0 * PI) * self.length_scale;
        [
            scalar::zero::<T>(),
            two_pi_l, // x: [0, 2πL]
            scalar::zero::<T>(),
            two_pi_l, // y: [0, 2πL]
            scalar::zero::<T>(),
            if self.is_3d {
                two_pi_l
            } else {
                scalar::zero::<T>()
            }, // z
        ]
    }

    fn length_scale(&self) -> T {
        self.length_scale
    }

    fn velocity_scale(&self) -> T {
        self.velocity_scale
    }
}
