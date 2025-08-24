//! Analytical solutions for CFD validation.
//!
//! This module provides exact analytical solutions for various fluid flow problems
//! that can be used to validate numerical CFD solvers.

// Result type available from cfd_core when needed
use nalgebra::{RealField, Vector3};
use num_traits::cast::FromPrimitive;
use std::f64::consts::PI;

/// Trait for analytical solutions
pub trait AnalyticalSolution<T: RealField + Copy> {
    /// Evaluate the solution at given coordinates and time
    fn evaluate(&self, x: T, y: T, z: T, t: T) -> Vector3<T>;

    /// Get the velocity field at given coordinates and time
    fn velocity(&self, x: T, y: T, z: T, t: T) -> Vector3<T> {
        self.evaluate(x, y, z, t)
    }

    /// Get the pressure field at given coordinates and time
    fn pressure(&self, x: T, y: T, z: T, t: T) -> T;

    /// Get the name of the analytical solution
    fn name(&self) -> &str;

    /// Get the domain bounds [`x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max`]
    fn domain_bounds(&self) -> [T; 6];

    /// Check if the solution is valid at given coordinates
    fn is_valid_at(&self, x: T, y: T, z: T) -> bool {
        let bounds = self.domain_bounds();
        x >= bounds[0] && x <= bounds[1] &&
        y >= bounds[2] && y <= bounds[3] &&
        z >= bounds[4] && z <= bounds[5]
    }
}

/// Poiseuille flow in a channel (1D/2D)
///
/// **Literature References:**
/// - Poiseuille, J.L.M. (1846). "Recherches expérimentales sur le mouvement des liquides"
/// - White, F.M. (2011). "Fluid Mechanics", 7th Edition, McGraw-Hill
/// - Schlichting, H. & Gersten, K. (2017). "Boundary-Layer Theory", 9th Edition
///
/// **Mathematical Formulation:**
/// - 2D channel: u(y) = `u_max` * (1 - (y/h)²) where h is half-height
/// - Pipe flow: u(r) = `u_max` * (1 - (r/R)²) where R is radius
pub struct PoiseuilleFlow<T: RealField + Copy> {
    /// Maximum velocity at channel center
    pub u_max: T,
    /// Channel half-width (for 2D) or radius (for cylindrical)
    pub channel_width: T,
    /// Pressure gradient (dp/dx)
    pub pressure_gradient: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Channel length
    pub length: T,
    /// Flow type: true for 2D channel, false for cylindrical pipe
    pub is_2d_channel: bool,
}

impl<T: RealField + Copy + FromPrimitive + Copy> PoiseuilleFlow<T> {
    /// Create new Poiseuille flow solution
    pub fn new(
        u_max: T,
        channel_width: T,
        pressure_gradient: T,
        viscosity: T,
        length: T,
        is_2d_channel: bool,
    ) -> Self {
        Self {
            u_max,
            channel_width,
            pressure_gradient,
            viscosity,
            length,
            is_2d_channel,
        }
    }

    /// Create 2D channel Poiseuille flow
    pub fn channel_2d(u_max: T, half_width: T, length: T, viscosity: T) -> Self {
        // For 2D channel: u_max = -dp/dx * h^2 / (2*mu)
        // So dp/dx = -2*mu*u_max / h^2
        let pressure_gradient = -T::from_f64(2.0).unwrap_or_else(|| T::zero()) * viscosity * u_max /
                               (half_width * half_width);

        Self::new(u_max, half_width, pressure_gradient, viscosity, length, true)
    }

    /// Create cylindrical pipe Poiseuille flow
    pub fn pipe_cylindrical(u_max: T, radius: T, length: T, viscosity: T) -> Self {
        // For cylindrical pipe: u_max = -dp/dx * R^2 / (4*mu)
        // So dp/dx = -4*mu*u_max / R^2
        let pressure_gradient = -T::from_f64(4.0).unwrap_or_else(|| T::zero()) * viscosity * u_max /
                               (radius * radius);

        Self::new(u_max, radius, pressure_gradient, viscosity, length, false)
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> AnalyticalSolution<T> for PoiseuilleFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        if self.is_2d_channel {
            // Parallel plate channel flow (y from 0 to h)
            // Exact solution: u(y) = -(1/2μ) * (dp/dx) * y * (h - y)
            // For normalized form with u_max at center:
            // u(y) = 4 * u_max * (y/h) * (1 - y/h)
            let y_norm = y / self.channel_width;
            let u = if y_norm >= T::zero() && y_norm <= T::one() {
                let coeff = T::from_f64(cfd_core::constants::physics::fluid::CHANNEL_FLOW_COEFFICIENT)
                    .unwrap_or_else(|| T::zero());
                coeff * self.u_max * y_norm * (T::one() - y_norm)
            } else {
                T::zero()
            };
            Vector3::new(u, T::zero(), T::zero())
        } else {
            // Hagen-Poiseuille flow in cylindrical pipe
            // u(r) = u_max * (1 - (r/R)^2)
            let r = (y * y + _z * _z).sqrt();
            let r_normalized = r / self.channel_width;
            let u = if r_normalized <= T::one() {
                self.u_max * (T::one() - r_normalized * r_normalized)
            } else {
                T::zero()
            };
            Vector3::new(u, T::zero(), T::zero())
        }
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        // Linear pressure drop: p(x) = p0 + (dp/dx) * x
        // Assuming p0 = 0 at x = 0
        self.pressure_gradient * x
    }

    fn name(&self) -> &str {
        if self.is_2d_channel {
            "Poiseuille Flow (2D Channel)"
        } else {
            "Poiseuille Flow (Cylindrical Pipe)"
        }
    }

    fn domain_bounds(&self) -> [T; 6] {
        if self.is_2d_channel {
            [T::zero(), self.length,
             -self.channel_width, self.channel_width,
             T::zero(), T::zero()]
        } else {
            [T::zero(), self.length,
             -self.channel_width, self.channel_width,
             -self.channel_width, self.channel_width]
        }
    }
}

/// Couette flow between parallel plates
///
/// **Literature References:**
/// - Couette, M. (1890). "Études sur le frottement des liquides"
/// - White, F.M. (2011). "Fluid Mechanics", 7th Edition, McGraw-Hill
/// - Kundu, P.K. & Cohen, I.M. (2016). "Fluid Mechanics", 6th Edition
///
/// **Mathematical Formulation:**
/// u(y) = U * y/h + (dp/dx) * y * (h - y) / (2*μ)
/// where U is plate velocity, h is gap, dp/dx is pressure gradient
pub struct CouetteFlow<T: RealField + Copy> {
    /// Velocity of the moving plate
    pub plate_velocity: T,
    /// Gap between plates
    pub gap: T,
    /// Pressure gradient (optional)
    pub pressure_gradient: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Length of the plates
    pub length: T,
}

impl<T: RealField + Copy + FromPrimitive + Copy> CouetteFlow<T> {
    /// Create new Couette flow solution
    pub fn new(plate_velocity: T, gap: T, pressure_gradient: T, viscosity: T, length: T) -> Self {
        Self {
            plate_velocity,
            gap,
            pressure_gradient,
            viscosity,
            length,
        }
    }

    /// Create Couette flow without pressure gradient
    pub fn no_pressure_gradient(plate_velocity: T, gap: T, length: T) -> Self {
        Self::new(plate_velocity, gap, T::zero(), T::one(), length)
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> AnalyticalSolution<T> for CouetteFlow<T> {
    fn evaluate(&self, _x: T, y: T, _z: T, _t: T) -> Vector3<T> {
        // Couette flow: u(y) = U * y/h + (dp/dx) * y * (h - y) / (2*mu)
        let y_normalized = y / self.gap;
        let linear_term = self.plate_velocity * y_normalized;

        let pressure_term = if self.pressure_gradient != T::zero() {
            let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
            self.pressure_gradient * y * (self.gap - y) /
            (two * self.viscosity)
        } else {
            T::zero()
        };

        let u = linear_term + pressure_term;
        Vector3::new(u, T::zero(), T::zero())
    }

    fn pressure(&self, x: T, _y: T, _z: T, _t: T) -> T {
        self.pressure_gradient * x
    }

    fn name(&self) -> &str {
        "Couette Flow"
    }

    fn domain_bounds(&self) -> [T; 6] {
        [T::zero(), self.length,
         T::zero(), self.gap,
         T::zero(), T::zero()]
    }
}

/// Taylor-Green vortex solution
pub struct TaylorGreenVortex<T: RealField + Copy> {
    /// Amplitude of the vortex
    pub amplitude: T,
    /// Kinematic viscosity
    pub viscosity: T,
    /// Domain size (assumed square)
    pub domain_size: T,
}

impl<T: RealField + Copy + FromPrimitive + Copy> TaylorGreenVortex<T> {
    /// Create new Taylor-Green vortex solution
    pub fn new(amplitude: T, viscosity: T, domain_size: T) -> Self {
        Self {
            amplitude,
            viscosity,
            domain_size,
        }
    }

    /// Calculate kinetic energy at time t
    pub fn kinetic_energy(&self, t: T) -> T {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        
        // Decay factor: exp(-4*nu*t) for kinetic energy
        let decay = (-T::from_f64(4.0).unwrap_or_else(|| T::zero()) * self.viscosity * t).exp();
        
        // Kinetic energy = 0.5 * amplitude^2 * exp(-4*nu*t)
        half * self.amplitude * self.amplitude * decay
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> AnalyticalSolution<T> for TaylorGreenVortex<T> {
    fn evaluate(&self, x: T, y: T, _z: T, t: T) -> Vector3<T> {
        let pi = T::from_f64(PI).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Decay factor: exp(-2*nu*t)
        let decay = (-two * self.viscosity * t).exp();

        // Normalized coordinates
        let kx = two * pi * x / self.domain_size;
        let ky = two * pi * y / self.domain_size;

        // Velocity components
        let u = self.amplitude * kx.cos() * ky.sin() * decay;
        let v = -self.amplitude * kx.sin() * ky.cos() * decay;

        Vector3::new(u, v, T::zero())
    }

    fn pressure(&self, x: T, y: T, _z: T, t: T) -> T {
        let pi = T::from_f64(PI).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());

        // Decay factor: exp(-4*nu*t)
        let decay = (-four * self.viscosity * t).exp();

        // Normalized coordinates
        let kx = two * pi * x / self.domain_size;
        let ky = two * pi * y / self.domain_size;

        // Pressure field
        let amp_squared = self.amplitude * self.amplitude;
        -amp_squared * (kx.cos() * two + ky.cos() * two) * decay / four
    }

    fn name(&self) -> &str {
        "Taylor-Green Vortex"
    }

    fn domain_bounds(&self) -> [T; 6] {
        [T::zero(), self.domain_size,
         T::zero(), self.domain_size,
         T::zero(), T::zero()]
    }
}

/// Stokes flow around a sphere
pub struct StokesFlow<T: RealField + Copy> {
    /// Sphere radius
    pub radius: T,
    /// Free stream velocity
    pub u_infinity: T,
    /// Dynamic viscosity
    pub viscosity: T,
    /// Sphere center coordinates
    pub center: Vector3<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> StokesFlow<T> {
    /// Create new Stokes flow solution
    pub fn new(radius: T, u_infinity: T, viscosity: T, center: Vector3<T>) -> Self {
        Self {
            radius,
            u_infinity,
            viscosity,
            center,
        }
    }

    /// Create Stokes flow with sphere at origin
    pub fn at_origin(radius: T, u_infinity: T, viscosity: T) -> Self {
        Self::new(radius, u_infinity, viscosity, Vector3::zeros())
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> AnalyticalSolution<T> for StokesFlow<T> {
    fn evaluate(&self, x: T, y: T, z: T, _t: T) -> Vector3<T> {
        // Position relative to sphere center
        let pos = Vector3::new(x, y, z) - self.center;
        let r = pos.norm();

        if r <= self.radius {
            // Inside sphere: zero velocity
            return Vector3::zeros();
        }

        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());

        // Stokes flow solution
        let a_over_r = self.radius / r;
        let a_over_r_cubed = a_over_r * a_over_r * a_over_r;

        // Velocity components (assuming flow in x-direction)
        let cos_theta = pos.x / r;
        let sin_theta = (pos.y * pos.y + pos.z * pos.z).sqrt() / r;

        let u_r = self.u_infinity * cos_theta *
                  (T::one() - three * a_over_r / T::from_f64(2.0).unwrap_or_else(|| T::zero()) + a_over_r_cubed / T::from_f64(2.0).unwrap_or_else(|| T::zero()));
        let u_theta = -self.u_infinity * sin_theta *
                      (T::one() - three * a_over_r / four - a_over_r_cubed / four);

        // Convert to Cartesian coordinates
        let u_x = u_r * cos_theta - u_theta * sin_theta;
        let u_y = if sin_theta != T::zero() {
            (u_r * pos.y / r + u_theta * pos.y / (r * sin_theta)) / r
        } else {
            T::zero()
        };
        let u_z = if sin_theta != T::zero() {
            (u_r * pos.z / r + u_theta * pos.z / (r * sin_theta)) / r
        } else {
            T::zero()
        };

        Vector3::new(u_x, u_y, u_z)
    }

    fn pressure(&self, x: T, y: T, z: T, _t: T) -> T {
        // Position relative to sphere center
        let pos = Vector3::new(x, y, z) - self.center;
        let r = pos.norm();

        if r <= self.radius {
            // Inside sphere: constant pressure
            return T::zero();
        }

        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Pressure field
        let cos_theta = pos.x / r;
        -three * self.viscosity * self.u_infinity * self.radius * cos_theta /
         (two * r * r)
    }

    fn name(&self) -> &str {
        "Stokes Flow Around Sphere"
    }

    fn domain_bounds(&self) -> [T; 6] {
        let bound = T::from_f64(10.0).unwrap_or_else(|| T::zero()) * self.radius;
        [self.center.x - bound, self.center.x + bound,
         self.center.y - bound, self.center.y + bound,
         self.center.z - bound, self.center.z + bound]
    }
}

/// Utility functions for analytical solutions
pub struct AnalyticalUtils;

impl AnalyticalUtils {
    /// Create a grid of points for evaluation
    pub fn create_grid<T: RealField + Copy + FromPrimitive + Copy>(
        bounds: [T; 6],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Vec<(T, T, T)> {
        let mut points = Vec::with_capacity(nx * ny * nz);

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = if nx > 1 {
                        bounds[0] + (bounds[1] - bounds[0]) *
                        T::from_usize(i).unwrap_or_else(|| T::zero()) / T::from_usize(nx - 1).unwrap_or_else(|| T::zero())
                    } else {
                        (bounds[0] + bounds[1]) / T::from_f64(2.0).unwrap_or_else(|| T::zero())
                    };

                    let y = if ny > 1 {
                        bounds[2] + (bounds[3] - bounds[2]) *
                        T::from_usize(j).unwrap_or_else(|| T::zero()) / T::from_usize(ny - 1).unwrap_or_else(|| T::zero())
                    } else {
                        (bounds[2] + bounds[3]) / T::from_f64(2.0).unwrap_or_else(|| T::zero())
                    };

                    let z = if nz > 1 {
                        bounds[4] + (bounds[5] - bounds[4]) *
                        T::from_usize(k).unwrap_or_else(|| T::zero()) / T::from_usize(nz - 1).unwrap_or_else(|| T::zero())
                    } else {
                        (bounds[4] + bounds[5]) / T::from_f64(2.0).unwrap_or_else(|| T::zero())
                    };

                    points.push((x, y, z));
                }
            }
        }

        points
    }

    /// Evaluate solution on a grid
    pub fn evaluate_on_grid<T, S>(
        solution: &S,
        points: &[(T, T, T)],
        time: T,
    ) -> (Vec<Vector3<T>>, Vec<T>)
    where
        T: RealField + Copy,
        S: AnalyticalSolution<T>,
    {
        let velocities: Vec<Vector3<T>> = points
            .iter()
            .map(|(x, y, z)| solution.velocity(*x, *y, *z, time))
            .collect();

        let pressures: Vec<T> = points
            .iter()
            .map(|(x, y, z)| solution.pressure(*x, *y, *z, time))
            .collect();

        (velocities, pressures)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_poiseuille_flow_2d() {
        let flow = PoiseuilleFlow::channel_2d(1.0, 1.0, 10.0, 0.001);

        // Test centerline velocity (at y = h/2 = 0.5)
        let vel_center = flow.velocity(5.0, 0.5, 0.0, 0.0);
        assert_relative_eq!(vel_center.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(vel_center.y, 0.0, epsilon = 1e-10);

        // Test wall velocities (should be zero at y=0 and y=1)
        let vel_wall_bottom = flow.velocity(5.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel_wall_bottom.x, 0.0, epsilon = 1e-10);
        let vel_wall_top = flow.velocity(5.0, 1.0, 0.0, 0.0);
        assert_relative_eq!(vel_wall_top.x, 0.0, epsilon = 1e-10);

        // Test pressure gradient
        let p1 = flow.pressure(0.0, 0.0, 0.0, 0.0);
        let p2 = flow.pressure(1.0, 0.0, 0.0, 0.0);
        let dp_dx = p2 - p1;
        assert!(dp_dx < 0.0); // Pressure should decrease downstream

        // Test domain bounds
        let bounds = flow.domain_bounds();
        assert_eq!(bounds[0], 0.0); // x_min
        assert_eq!(bounds[1], 10.0); // x_max
        assert_eq!(bounds[2], -1.0); // y_min
        assert_eq!(bounds[3], 1.0); // y_max
    }

    #[test]
    fn test_poiseuille_flow_cylindrical() {
        let flow = PoiseuilleFlow::pipe_cylindrical(2.0, 0.5, 5.0, 0.001);

        // Test centerline velocity
        let vel_center = flow.velocity(2.5, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel_center.x, 2.0, epsilon = 1e-10);

        // Test velocity at r = R/2
        let vel_half = flow.velocity(2.5, 0.25, 0.0, 0.0);
        let expected = 2.0 * (1.0 - 0.25); // u_max * (1 - (r/R)^2)
        assert_relative_eq!(vel_half.x, expected, epsilon = 1e-10);

        // Test velocity at wall (r = R)
        let vel_wall = flow.velocity(2.5, 0.5, 0.0, 0.0);
        assert_relative_eq!(vel_wall.x, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_couette_flow() {
        let flow = CouetteFlow::no_pressure_gradient(1.0, 2.0, 10.0);

        // Test velocity at bottom plate
        let vel_bottom = flow.velocity(5.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel_bottom.x, 0.0, epsilon = 1e-10);

        // Test velocity at top plate
        let vel_top = flow.velocity(5.0, 2.0, 0.0, 0.0);
        assert_relative_eq!(vel_top.x, 1.0, epsilon = 1e-10);

        // Test velocity at middle
        let vel_middle = flow.velocity(5.0, 1.0, 0.0, 0.0);
        assert_relative_eq!(vel_middle.x, 0.5, epsilon = 1e-10);
    }

    #[test]
    fn test_taylor_green_vortex() {
        let vortex = TaylorGreenVortex::new(1.0, 0.01, 2.0 * PI);

        // Test velocity at t=0, x=0, y=0 (should be zero due to sin(0) and cos(0))
        let vel = vortex.velocity(0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel.x, 0.0, epsilon = 1e-10); // cos(0) * sin(0) = 0
        assert_relative_eq!(vel.y, 0.0, epsilon = 1e-10); // -sin(0) * cos(0) = 0

        // Test pressure at origin
        let p = vortex.pressure(0.0, 0.0, 0.0, 0.0);
        assert!(p < 0.0); // Should be negative due to vortex

        // Test decay over time
        let vel_t0 = vortex.velocity(PI / 2.0, PI / 2.0, 0.0, 0.0);
        let vel_t1 = vortex.velocity(PI / 2.0, PI / 2.0, 0.0, 1.0);
        assert!(vel_t1.norm() < vel_t0.norm()); // Should decay over time
    }

    #[test]
    fn test_stokes_flow() {
        let flow = StokesFlow::at_origin(1.0, 1.0, 0.001);

        // Test velocity inside sphere (should be zero)
        let vel_inside = flow.velocity(0.5, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel_inside.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(vel_inside.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(vel_inside.z, 0.0, epsilon = 1e-10);

        // Test velocity far from sphere (should approach free stream)
        let vel_far = flow.velocity(100.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(vel_far.x, 1.0, epsilon = 2e-2); // Should approach u_infinity (looser tolerance)

        // Test that velocity is finite at sphere surface
        let vel_surface = flow.velocity(1.0_f64, 0.0, 0.0, 0.0);
        assert!(vel_surface.x.is_finite());
        assert!(vel_surface.y.is_finite());
        assert!(vel_surface.z.is_finite());
    }

    #[test]
    fn test_analytical_utils_grid() {
        let bounds = [0.0, 1.0, 0.0, 1.0, 0.0, 0.0];
        let points = AnalyticalUtils::create_grid(bounds, 3, 3, 1);

        assert_eq!(points.len(), 9);

        // Check corner points
        assert_relative_eq!(points[0].0, 0.0, epsilon = 1e-10);
        assert_relative_eq!(points[0].1, 0.0, epsilon = 1e-10);
        assert_relative_eq!(points[8].0, 1.0, epsilon = 1e-10);
        assert_relative_eq!(points[8].1, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_evaluate_on_grid() {
        let flow = PoiseuilleFlow::channel_2d(1.0, 1.0, 2.0, 0.001);
        let bounds = flow.domain_bounds();
        let points = AnalyticalUtils::create_grid(bounds, 3, 3, 1);

        let (velocities, pressures) = AnalyticalUtils::evaluate_on_grid(&flow, &points, 0.0);

        assert_eq!(velocities.len(), 9);
        assert_eq!(pressures.len(), 9);

        // All velocities should be in x-direction only
        for vel in &velocities {
            assert_relative_eq!(vel.y, 0.0, epsilon = 1e-10);
            assert_relative_eq!(vel.z, 0.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_domain_validation() {
        let flow = PoiseuilleFlow::channel_2d(1.0, 1.0, 10.0, 0.001);

        // Test valid points
        assert!(flow.is_valid_at(5.0, 0.0, 0.0));
        assert!(flow.is_valid_at(0.0, -1.0, 0.0));
        assert!(flow.is_valid_at(10.0, 1.0, 0.0));

        // Test invalid points
        assert!(!flow.is_valid_at(-1.0, 0.0, 0.0)); // x < 0
        assert!(!flow.is_valid_at(11.0, 0.0, 0.0)); // x > length
        assert!(!flow.is_valid_at(5.0, 2.0, 0.0)); // y > channel_width
        assert!(!flow.is_valid_at(5.0, -2.0, 0.0)); // y < -channel_width
    }
}
