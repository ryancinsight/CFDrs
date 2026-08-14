//! Stokes flow solutions - low Reynolds number flows

use super::AnalyticalSolution;
use crate::scalar;
use aequitas::systems::si::quantities::{
    Angle, Dimensionless, DynamicViscosity, Force, Length, MassDensity, Velocity,
    VolumetricFlowRate,
};
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector3;
use std::f64::consts::PI;

/// Stokes flow around a sphere
///
/// Analytical solution for creeping flow (Re << 1) around a sphere.
pub struct StokesFlow<T: RealField + Copy> {
    /// Sphere radius
    pub sphere_radius: Length<T>,
    /// Free stream velocity
    pub free_stream_velocity: Velocity<T>,
    /// Dynamic viscosity
    pub viscosity: DynamicViscosity<T>,
    /// Fluid density
    pub density: MassDensity<T>,
}

impl<T: RealField + Copy + FloatElement> StokesFlow<T> {
    /// Create Stokes flow solution
    pub fn create(
        sphere_radius: Length<T>,
        free_stream_velocity: Velocity<T>,
        viscosity: DynamicViscosity<T>,
        density: MassDensity<T>,
    ) -> Self {
        Self {
            sphere_radius,
            free_stream_velocity,
            viscosity,
            density,
        }
    }

    /// Get the drag force on the sphere (Stokes' law)
    pub fn drag_force(&self) -> Force<T> {
        let pi = <T as FloatElement>::from_f64(PI);
        let six = <T as FloatElement>::from_f64(6.0);

        Force::from_base(
            six * pi
                * self.viscosity.into_base()
                * self.sphere_radius.into_base()
                * self.free_stream_velocity.into_base(),
        )
    }

    /// Get the drag coefficient
    pub fn drag_coefficient(&self) -> Dimensionless<T> {
        let reynolds = self.reynolds_number().into_base();
        if reynolds > <T as FloatElement>::from_f64(0.01) {
            // Stokes drag coefficient: CD = 24/Re
            Dimensionless::from_base(<T as FloatElement>::from_f64(24.0) / reynolds)
        } else {
            // Avoid division by very small number
            Dimensionless::from_base(<T as FloatElement>::from_f64(2400.0))
        }
    }

    /// Get Reynolds number based on sphere diameter
    pub fn reynolds_number(&self) -> Dimensionless<T> {
        let diameter = <T as FloatElement>::from_f64(2.0) * self.sphere_radius.into_base();
        Dimensionless::from_base(
            self.density.into_base() * self.free_stream_velocity.into_base() * diameter
                / self.viscosity.into_base(),
        )
    }

    /// Get the stream function value at (r, θ)
    pub fn stream_function(&self, r: Length<T>, theta: Angle<T>) -> VolumetricFlowRate<T> {
        let r = r.into_base();
        let a = self.sphere_radius.into_base();
        let u_inf = self.free_stream_velocity.into_base();

        // ψ = (U∞/2) * r² * sin²(θ) * (1 - 3a/(2r) + a³/(2r³))
        let sin_theta = scalar::sin(theta.into_base());
        let sin2_theta = sin_theta * sin_theta;

        let half = <T as FloatElement>::from_f64(0.5);
        let three_halves = <T as FloatElement>::from_f64(1.5);

        let term1 = scalar::one::<T>();
        let term2 = -three_halves * a / r;
        let term3 = half * a * a * a / (r * r * r);

        VolumetricFlowRate::from_base(half * u_inf * r * r * sin2_theta * (term1 + term2 + term3))
    }
}

impl<T: RealField + Copy + FloatElement> AnalyticalSolution<T> for StokesFlow<T> {
    fn evaluate(&self, x: T, y: T, z: T, _t: T) -> Vector3<T> {
        // Convert to spherical coordinates (r, θ, φ)
        let r = scalar::sqrt(x * x + y * y + z * z);

        let sphere_radius = self.sphere_radius.into_base();
        let free_stream_velocity = self.free_stream_velocity.into_base();

        if r < sphere_radius {
            // Inside sphere: no flow
            return Vector3::zeros();
        }

        let a = sphere_radius;
        let u_inf = free_stream_velocity;

        // For simplicity, assume flow in x-direction
        // u_r = U∞ * cos(θ) * (1 - 3a/(2r) + a³/(2r³))
        // u_θ = -U∞ * sin(θ) * (1 - 3a/(4r) - a³/(4r³))

        let cos_theta = if r > scalar::zero::<T>() {
            x / r
        } else {
            scalar::zero::<T>()
        };
        let sin_theta = if r > scalar::zero::<T>() {
            scalar::sqrt(y * y + z * z) / r
        } else {
            scalar::zero::<T>()
        };

        let three_halves = <T as FloatElement>::from_f64(1.5);
        let half = <T as FloatElement>::from_f64(0.5);
        let quarter = <T as FloatElement>::from_f64(0.25);

        // Radial velocity component
        let u_r = u_inf
            * cos_theta
            * (scalar::one::<T>() - three_halves * a / r + half * a * a * a / (r * r * r));

        // Tangential velocity component
        let u_theta = -u_inf
            * sin_theta
            * (scalar::one::<T>()
                - three_halves * quarter * a / r
                - quarter * a * a * a / (r * r * r));

        // Convert from spherical (r, θ, φ) to Cartesian coordinates (x, y, z)
        // Full transformation with proper handling of all spherical angles
        let phi = scalar::atan2(z, y); // Azimuthal angle in x-y plane
        let cos_phi = scalar::cos(phi);
        let sin_phi = scalar::sin(phi);

        // Complete spherical to Cartesian velocity transformation
        let u_x = u_r * sin_theta * cos_phi - u_theta * cos_theta * cos_phi;
        let u_y = u_r * sin_theta * sin_phi - u_theta * cos_theta * sin_phi;
        let u_z = u_r * cos_theta + u_theta * sin_theta;

        Vector3::new(u_x, u_y, u_z)
    }

    fn pressure(&self, x: T, y: T, z: T, _t: T) -> T {
        let r = scalar::sqrt(x * x + y * y + z * z);

        if r < self.sphere_radius.into_base() {
            // Inside sphere
            return scalar::zero::<T>();
        }

        // p = p∞ - (3μU∞a/2) * (x/r³)
        let three_halves = <T as FloatElement>::from_f64(1.5);
        let pressure_drop = three_halves
            * self.viscosity.into_base()
            * self.free_stream_velocity.into_base()
            * self.sphere_radius.into_base()
            * x
            / (r * r * r);

        -pressure_drop // Relative to p∞
    }

    fn name(&self) -> &'static str {
        "Stokes Flow Around Sphere"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let domain_size = <T as FloatElement>::from_f64(10.0) * self.sphere_radius.into_base();
        [
            -domain_size,
            domain_size, // x
            -domain_size,
            domain_size, // y
            -domain_size,
            domain_size, // z
        ]
    }

    fn length_scale(&self) -> T {
        self.sphere_radius.into_base()
    }

    fn velocity_scale(&self) -> T {
        self.free_stream_velocity.into_base()
    }
}
