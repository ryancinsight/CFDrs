//! Manufactured solutions for the advection equation
//!
//! Solves: ∂u/∂t + v·∇u = S(x,y,t)
//! where S is the manufactured source term

use super::ManufacturedSolution;
use crate::scalar;
use eunomia::{FloatElement, RealField};
use std::f64::consts::PI;

/// Manufactured solution for 2D linear advection
#[derive(Clone, Copy)]
pub struct ManufacturedAdvection<T: RealField + Copy> {
    /// Velocity in x-direction
    pub vx: T,
    /// Velocity in y-direction
    pub vy: T,
}

impl<T: RealField + Copy + FloatElement> ManufacturedAdvection<T> {
    /// Create a new manufactured advection solution
    pub fn new(vx: T, vy: T) -> Self {
        Self { vx, vy }
    }

    /// Create with unit velocity
    pub fn unit_velocity() -> Self {
        Self {
            vx: scalar::one::<T>(),
            vy: scalar::one::<T>(),
        }
    }
}

impl<T: RealField + Copy + FloatElement> ManufacturedSolution<T> for ManufacturedAdvection<T> {
    /// Exact solution: u(x,y,t) = sin(x - vx*t) * cos(y - vy*t)
    /// This represents a wave propagating with velocity (vx, vy)
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        scalar::sin(x - self.vx * t) * scalar::cos(y - self.vy * t)
    }

    /// Source term: S = ∂u/∂t + v·∇u
    /// For pure advection with no source, this should be zero
    /// We add a manufactured source to test the solver
    fn source_term(&self, _x: T, _y: T, _z: T, _t: T) -> T {
        // For pure advection, the source term is zero
        // But we can add a manufactured source for testing
        scalar::zero::<T>()
    }
}

/// Manufactured solution for rotating advection (solid body rotation)
pub struct RotatingAdvection<T: RealField + Copy> {
    /// Angular velocity
    pub omega: T,
    /// Center of rotation x-coordinate
    pub cx: T,
    /// Center of rotation y-coordinate
    pub cy: T,
}

impl<T: RealField + Copy + FloatElement> RotatingAdvection<T> {
    /// Create a new rotating advection solution
    pub fn new(omega: T) -> Self {
        Self {
            omega,
            cx: scalar::from_f64(0.5),
            cy: scalar::from_f64(0.5),
        }
    }

    /// Get velocity at a point for solid body rotation
    pub fn velocity(&self, x: T, y: T) -> (T, T) {
        let dx = x - self.cx;
        let dy = y - self.cy;
        let vx = -self.omega * dy;
        let vy = self.omega * dx;
        (vx, vy)
    }
}

impl<T: RealField + Copy + FloatElement> ManufacturedSolution<T> for RotatingAdvection<T> {
    /// Exact solution: Gaussian bump rotating around center
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // Rotate the coordinates backward in time
        let theta = -self.omega * t;
        let cos_theta = scalar::cos(theta);
        let sin_theta = scalar::sin(theta);

        // Transform to rotated coordinates
        let dx = x - self.cx;
        let dy = y - self.cy;
        let x_rot = dx * cos_theta - dy * sin_theta;
        let y_rot = dx * sin_theta + dy * cos_theta;

        // Gaussian bump at rotated position
        let sigma = scalar::from_f64::<T>(0.1);
        let r_squared = x_rot * x_rot + y_rot * y_rot;
        scalar::exp(-r_squared / (scalar::from_f64::<T>(2.0) * sigma * sigma))
    }

    /// Source term for rotating advection
    fn source_term(&self, _x: T, _y: T, _z: T, _t: T) -> T {
        // For solid body rotation with no diffusion, source is zero
        scalar::zero::<T>()
    }
}

/// Manufactured solution for Burgers' equation
pub struct ManufacturedBurgers<T: RealField + Copy> {
    /// Viscosity coefficient
    pub nu: T,
}

impl<T: RealField + Copy + FloatElement> ManufacturedBurgers<T> {
    /// Create a new manufactured Burgers solution
    pub fn new(nu: T) -> Self {
        Self { nu }
    }
}

impl<T: RealField + Copy + FloatElement> ManufacturedSolution<T> for ManufacturedBurgers<T> {
    /// Exact solution: u(x,t) = 2*nu*π*sin(π*x)*exp(-nu*π²*t) / (2 + cos(π*x)*exp(-nu*π²*t))
    /// This is the Cole-Hopf solution for Burgers' equation
    fn exact_solution(&self, x: T, _y: T, _z: T, t: T) -> T {
        let pi = scalar::from_f64::<T>(PI);
        let two = scalar::from_f64::<T>(2.0);
        let exp_term = scalar::exp(-self.nu * pi * pi * t);
        let numerator = two * self.nu * pi * scalar::sin(pi * x) * exp_term;
        let denominator = two + scalar::cos(pi * x) * exp_term;
        numerator / denominator
    }

    /// Source term for Burgers' equation: ∂u/∂t + u*∂u/∂x = nu*∂²u/∂x²
    fn source_term(&self, _x: T, _y: T, _z: T, _t: T) -> T {
        // For the Cole-Hopf solution, no source term is needed
        scalar::zero::<T>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_advection_solution() {
        let solution = ManufacturedAdvection::new(1.0_f64, 1.0);

        // Test wave propagation
        let u0 = solution.exact_solution(0.0, 0.0, 0.0, 0.0);
        let u1 = solution.exact_solution(1.0, 1.0, 0.0, 1.0);
        assert!((u0 - u1).abs() < 1e-10);
    }

    #[test]
    fn test_rotating_advection() {
        let solution = RotatingAdvection::new(2.0_f64 * PI);

        // After one full rotation (t=1), solution should return to initial state
        let u0 = solution.exact_solution(0.6, 0.5, 0.0, 0.0);
        let u1 = solution.exact_solution(0.6, 0.5, 0.0, 1.0);
        assert!((u0 - u1).abs() < 1e-10);
    }

    #[test]
    fn test_burgers_decay() {
        let solution = ManufacturedBurgers::new(0.01_f64);

        // Solution should decay with time
        let u0 = solution.exact_solution(0.5, 0.0, 0.0, 0.0);
        let u1 = solution.exact_solution(0.5, 0.0, 0.0, 1.0);
        assert!(u1.abs() < u0.abs());
    }
}
