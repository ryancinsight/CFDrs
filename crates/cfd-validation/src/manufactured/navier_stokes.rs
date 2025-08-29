//! Manufactured solutions for the Navier-Stokes equations
//!
//! Solves the incompressible Navier-Stokes equations:
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f
//! ∇·u = 0

use super::ManufacturedSolution;
use nalgebra::{RealField, Vector2};
use num_traits::Float;
use std::f64::consts::PI;

/// Manufactured solution for 2D incompressible Navier-Stokes
pub struct ManufacturedNavierStokes<T: RealField + Float> {
    /// Kinematic viscosity
    pub nu: T,
    /// Density
    pub rho: T,
    /// Characteristic length
    pub length: T,
    /// Characteristic velocity
    pub velocity: T,
}

impl<T: RealField + Float> ManufacturedNavierStokes<T> {
    /// Create a new manufactured Navier-Stokes solution
    pub fn new(nu: T, rho: T) -> Self {
        Self {
            nu,
            rho,
            length: T::one(),
            velocity: T::one(),
        }
    }

    /// Get the Reynolds number
    pub fn reynolds_number(&self) -> T {
        self.velocity * self.length / self.nu
    }
}

impl<T: RealField + Float> ManufacturedSolution<T> for ManufacturedNavierStokes<T> {
    /// Exact solution for velocity and pressure
    /// Using the Taylor-Green vortex as manufactured solution
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // This returns the u-component of velocity
        // For full solution, we'd need separate methods for u, v, and p
        let pi = T::from(PI).unwrap();
        let decay = Float::exp(-T::from(2.0).unwrap() * self.nu * pi * pi * t);
        Float::sin(pi * x) * Float::cos(pi * y) * decay
    }

    /// Source term for the momentum equation
    fn source_term(&self, _x: T, _y: T, _z: T, _t: T) -> T {
        // For Taylor-Green vortex, the source term is zero
        // as it's an exact solution to the Navier-Stokes equations
        T::zero()
    }
}

/// Taylor-Green vortex solution for 2D Navier-Stokes
pub struct TaylorGreenManufactured<T: RealField + Float> {
    /// Kinematic viscosity
    pub nu: T,
    /// Wave number
    pub k: T,
}

impl<T: RealField + Float> TaylorGreenManufactured<T> {
    /// Create a new Taylor-Green manufactured solution
    pub fn new(nu: T) -> Self {
        let pi = T::from(PI).unwrap();
        Self { nu, k: pi }
    }

    /// Get velocity components
    pub fn velocity(&self, x: T, y: T, t: T) -> Vector2<T> {
        let decay = Float::exp(-T::from(2.0).unwrap() * self.nu * self.k * self.k * t);
        let u = Float::sin(self.k * x) * Float::cos(self.k * y) * decay;
        let v = -Float::cos(self.k * x) * Float::sin(self.k * y) * decay;
        Vector2::new(u, v)
    }

    /// Get pressure field
    pub fn pressure(&self, x: T, y: T, t: T) -> T {
        let decay = Float::exp(-T::from(4.0).unwrap() * self.nu * self.k * self.k * t);
        let quarter = T::from(0.25).unwrap();
        -quarter
            * (Float::cos(T::from(2.0).unwrap() * self.k * x)
                + Float::cos(T::from(2.0).unwrap() * self.k * y))
            * decay
    }

    /// Get vorticity
    pub fn vorticity(&self, x: T, y: T, t: T) -> T {
        let decay = Float::exp(-T::from(2.0).unwrap() * self.nu * self.k * self.k * t);
        -T::from(2.0).unwrap() * self.k * Float::sin(self.k * x) * Float::sin(self.k * y) * decay
    }

    /// Get kinetic energy
    pub fn kinetic_energy(&self, t: T) -> T {
        let decay = Float::exp(-T::from(4.0).unwrap() * self.nu * self.k * self.k * t);
        T::from(0.25).unwrap() * decay
    }

    /// Get enstrophy (integral of vorticity squared)
    pub fn enstrophy(&self, t: T) -> T {
        let decay = Float::exp(-T::from(4.0).unwrap() * self.nu * self.k * self.k * t);
        self.k * self.k * decay
    }
}

/// Kovasznay flow - exact solution for 2D steady Navier-Stokes
pub struct KovasznayFlow<T: RealField + Float> {
    /// Reynolds number
    pub re: T,
    /// Lambda parameter
    pub lambda: T,
}

impl<T: RealField + Float> KovasznayFlow<T> {
    /// Create a new Kovasznay flow solution
    pub fn new(re: T) -> Self {
        let half_re = re / T::from(2.0).unwrap();
        let pi_val = T::from(PI).unwrap();
        let discriminant = Float::sqrt(half_re * half_re + T::from(4.0).unwrap() * pi_val * pi_val);
        let lambda = half_re - discriminant;
        Self { re, lambda }
    }

    /// Get velocity components
    pub fn velocity(&self, x: T, y: T) -> Vector2<T> {
        let pi = T::from(PI).unwrap();
        let two_pi = T::from(2.0).unwrap() * pi;

        let exp_lambda_x = Float::exp(self.lambda * x);
        let u = T::one() - exp_lambda_x * Float::cos(two_pi * y);
        let v = self.lambda / two_pi * exp_lambda_x * Float::sin(two_pi * y);

        Vector2::new(u, v)
    }

    /// Get pressure field
    pub fn pressure(&self, x: T) -> T {
        let half = T::from(0.5).unwrap();
        half * (T::one() - Float::exp(T::from(2.0).unwrap() * self.lambda * x))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_taylor_green_decay() {
        let solution = TaylorGreenManufactured::new(0.01);

        // Kinetic energy should decay exponentially
        let ke0 = solution.kinetic_energy(0.0);
        let ke1 = solution.kinetic_energy(1.0);
        assert!(ke1 < ke0);

        // Check decay rate
        let expected_ratio = Float::exp(-4.0 * 0.01 * PI * PI);
        assert!((ke1 / ke0 - expected_ratio).abs() < 1e-10);
    }

    #[test]
    fn test_kovasznay_flow() {
        let solution = KovasznayFlow::new(40.0);

        // At x=0, y=0, velocity should be (1, 0)
        let vel = solution.velocity(0.0, 0.0);
        assert!((vel.x - 0.0).abs() < 1e-10);
        assert!(vel.y.abs() < 1e-10);
    }

    #[test]
    fn test_incompressibility() {
        let solution = TaylorGreenManufactured::new(0.01);

        // Check divergence-free condition
        let h = 0.001;
        let x = 0.5;
        let y = 0.5;
        let t = 0.0;

        let v_center = solution.velocity(x, y, t);
        let v_right = solution.velocity(x + h, y, t);
        let v_top = solution.velocity(x, y + h, t);

        let dudx = (v_right.x - v_center.x) / h;
        let dvdy = (v_top.y - v_center.y) / h;
        let divergence = dudx + dvdy;

        assert!(divergence.abs() < 1e-6);
    }
}
