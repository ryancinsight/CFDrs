//! Taylor-Green vortex - decaying vortex solution

use super::AnalyticalSolution;
use crate::scalar;
use aequitas::systems::si::quantities::{
    Dimensionless, Energy, Force, KinematicViscosity, Length, MassDensity, ReciprocalTime,
    ReciprocalTimeSquared, Time, Velocity,
};
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;
use std::f64::consts::PI;

/// Spatial dimensionality of the Taylor-Green solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TaylorGreenDimension {
    /// Two-dimensional flow with energy reported per unit depth.
    TwoD,
    /// Three-dimensional flow with volumetric energy.
    ThreeD,
}

/// Kinetic-energy result whose unit follows the solution dimensionality.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum TaylorGreenKineticEnergy<T> {
    /// Two-dimensional kinetic energy per unit depth (`N`).
    PerDepth(Force<T>),
    /// Three-dimensional kinetic energy (`J`).
    Volumetric(Energy<T>),
}

/// Taylor-Green vortex analytical solution.
///
/// Represents a decaying vortex flow that is an exact solution to the
/// incompressible Navier-Stokes equations in 2D/3D periodic domains.
pub struct TaylorGreenVortex<T: RealField + Copy> {
    /// Characteristic length scale.
    pub length_scale: Length<T>,
    /// Characteristic velocity scale.
    pub velocity_scale: Velocity<T>,
    /// Kinematic viscosity.
    pub viscosity: KinematicViscosity<T>,
    /// Density for pressure and energy calculations.
    pub density: MassDensity<T>,
    /// Spatial dimensionality of the solution.
    pub dimension: TaylorGreenDimension,
}

impl<T: RealField + Copy + FloatElement> TaylorGreenVortex<T> {
    /// Create a Taylor-Green vortex solution.
    pub fn create(
        length_scale: Length<T>,
        velocity_scale: Velocity<T>,
        viscosity: KinematicViscosity<T>,
        density: MassDensity<T>,
        dimension: TaylorGreenDimension,
    ) -> Self {
        Self {
            length_scale,
            velocity_scale,
            viscosity,
            density,
            dimension,
        }
    }

    /// Create a two-dimensional Taylor-Green vortex.
    pub fn create_2d(
        length_scale: Length<T>,
        velocity_scale: Velocity<T>,
        viscosity: KinematicViscosity<T>,
    ) -> Self {
        Self::create(
            length_scale,
            velocity_scale,
            viscosity,
            MassDensity::from_base(scalar::one::<T>()),
            TaylorGreenDimension::TwoD,
        )
    }

    /// Create a three-dimensional Taylor-Green vortex.
    pub fn create_3d(
        length_scale: Length<T>,
        velocity_scale: Velocity<T>,
        viscosity: KinematicViscosity<T>,
    ) -> Self {
        Self::create(
            length_scale,
            velocity_scale,
            viscosity,
            MassDensity::from_base(scalar::one::<T>()),
            TaylorGreenDimension::ThreeD,
        )
    }

    fn is_3d(&self) -> bool {
        self.dimension == TaylorGreenDimension::ThreeD
    }

    /// Get the Reynolds number.
    pub fn reynolds_number(&self) -> Dimensionless<T> {
        let value = self.velocity_scale.into_base() * self.length_scale.into_base()
            / self.viscosity.into_base();
        Dimensionless::from_base(value)
    }

    /// Get the decay rate.
    pub fn decay_rate(&self) -> ReciprocalTime<T> {
        let pi = scalar::from_f64::<T>(PI);
        let factor = if self.is_3d() {
            scalar::from_f64::<T>(3.0)
        } else {
            scalar::from_f64::<T>(2.0)
        };
        let length = self.length_scale.into_base();
        let value = factor * self.viscosity.into_base() * pi * pi / (length * length);
        ReciprocalTime::from_base(value)
    }

    /// Get kinetic energy at time `t`.
    pub fn kinetic_energy(&self, t: Time<T>) -> TaylorGreenKineticEnergy<T> {
        let length = self.length_scale.into_base();
        let velocity = self.velocity_scale.into_base();
        let density = self.density.into_base();
        let initial_energy = if self.is_3d() {
            // E₀ = (1/16) * ρ * U² * L³ for 3D.
            let factor = scalar::from_f64::<T>(1.0 / 16.0);
            TaylorGreenKineticEnergy::Volumetric(Energy::from_base(
                factor * density * velocity * velocity * length * length * length,
            ))
        } else {
            // E₀ = (1/4) * ρ * U² * L² for 2D, reported per unit depth.
            let factor = scalar::from_f64::<T>(0.25);
            TaylorGreenKineticEnergy::PerDepth(Force::from_base(
                factor * density * velocity * velocity * length * length,
            ))
        };

        let time = t.into_base();
        let decay =
            scalar::exp(-(scalar::from_f64::<T>(2.0) * self.decay_rate().into_base() * time));
        match initial_energy {
            TaylorGreenKineticEnergy::PerDepth(value) => {
                TaylorGreenKineticEnergy::PerDepth(Force::from_base(value.into_base() * decay))
            }
            TaylorGreenKineticEnergy::Volumetric(value) => {
                TaylorGreenKineticEnergy::Volumetric(Energy::from_base(value.into_base() * decay))
            }
        }
    }

    /// Get enstrophy (vorticity squared) at time `t`.
    pub fn enstrophy(&self, t: Time<T>) -> ReciprocalTimeSquared<T> {
        let pi = scalar::from_f64::<T>(PI);
        let velocity = self.velocity_scale.into_base();
        let length = self.length_scale.into_base();
        let initial_enstrophy = velocity * velocity * pi * pi / (length * length);
        let time = t.into_base();
        let decay =
            scalar::exp(-(scalar::from_f64::<T>(2.0) * self.decay_rate().into_base() * time));
        ReciprocalTimeSquared::from_base(initial_enstrophy * decay)
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for TaylorGreenVortex<T> {
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        let pi = scalar::from_f64::<T>(PI);
        let decay = scalar::exp(-self.decay_rate().into_base() * t);
        let length = self.length_scale.into_base();
        let velocity = self.velocity_scale.into_base();

        // Normalize coordinates at the mesh-coordinate boundary.
        let kx = pi * x / length;
        let ky = pi * y / length;

        if self.is_3d() {
            let kz = pi * z / length;

            // 3D Taylor-Green vortex.
            let u = velocity * scalar::sin(kx) * scalar::cos(ky) * scalar::cos(kz) * decay;
            let v = -velocity * scalar::cos(kx) * scalar::sin(ky) * scalar::cos(kz) * decay;
            let w = scalar::zero::<T>(); // Standard Taylor-Green case.

            Vector3::new(u, v, w)
        } else {
            // 2D Taylor-Green vortex.
            let u = velocity * scalar::cos(kx) * scalar::sin(ky) * decay;
            let v = -velocity * scalar::sin(kx) * scalar::cos(ky) * decay;

            Vector3::new(u, v, scalar::zero::<T>())
        }
    }

    fn pressure(&self, x: T, y: T, z: T, t: T) -> T {
        let pi = scalar::from_f64::<T>(PI);
        let decay = scalar::exp(-scalar::from_f64::<T>(2.0) * self.decay_rate().into_base() * t);
        let length = self.length_scale.into_base();
        let density = self.density.into_base();
        let velocity = self.velocity_scale.into_base();

        // Normalize coordinates at the mesh-coordinate boundary.
        let kx = pi * x / length;
        let ky = pi * y / length;

        if self.is_3d() {
            let kz = pi * z / length;
            let factor = scalar::from_f64::<T>(1.0 / 16.0);
            let two = scalar::from_f64::<T>(2.0);

            // p = ρU²/16 * (cos(2kx) + cos(2ky)) * (cos(2kz) + 2)
            //     * exp(-2νk²t).
            factor
                * density
                * velocity
                * velocity
                * (scalar::cos(two * kx) + scalar::cos(two * ky))
                * (scalar::cos(two * kz) + two)
                * decay
        } else {
            let factor = scalar::from_f64::<T>(0.25);
            let two = scalar::from_f64::<T>(2.0);

            // p = -ρU²/4 * (cos(2kx) + cos(2ky)) * exp(-2νk²t).
            -factor
                * density
                * velocity
                * velocity
                * (scalar::cos(two * kx) + scalar::cos(two * ky))
                * decay
        }
    }

    fn name(&self) -> &str {
        if self.is_3d() {
            "3D Taylor-Green Vortex"
        } else {
            "2D Taylor-Green Vortex"
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        let two_pi_l = scalar::from_f64::<T>(2.0 * PI) * self.length_scale.into_base();
        [
            scalar::zero::<T>(),
            two_pi_l, // x: [0, 2πL]
            scalar::zero::<T>(),
            two_pi_l, // y: [0, 2πL]
            scalar::zero::<T>(),
            if self.is_3d() {
                two_pi_l
            } else {
                scalar::zero::<T>()
            }, // z
        ]
    }

    fn length_scale(&self) -> T {
        self.length_scale.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.velocity_scale.into_base()
    }
}
