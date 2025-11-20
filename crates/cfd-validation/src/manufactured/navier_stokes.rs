//! Manufactured solutions for the Navier-Stokes equations
//!
//! Provides complete MMS (Method of Manufactured Solutions) for incompressible Navier-Stokes:
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f
//! ∇·u = 0
//!
//! References:
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"
//! - Salari, K. & Knupp, P. (2000) "Code Verification by the Method of Manufactured Solutions"

use super::ManufacturedSolution;
use cfd_core::conversion::SafeFromF64;
use nalgebra::{RealField, Vector2};
use num_traits::{Float, FromPrimitive};
use std::f64::consts::PI;

/// Complete manufactured solution for 2D incompressible Navier-Stokes equations
pub trait NavierStokesManufacturedSolution<T: RealField + Copy> {
    /// Exact velocity solution at (x,y,t)
    fn exact_velocity(&self, x: T, y: T, t: T) -> Vector2<T>;

    /// Exact pressure solution at (x,y,t)
    fn exact_pressure(&self, x: T, y: T, t: T) -> T;

    /// Source term for u-momentum equation
    fn momentum_source_u(&self, x: T, y: T, t: T) -> T;

    /// Source term for v-momentum equation
    fn momentum_source_v(&self, x: T, y: T, t: T) -> T;

    /// Verify continuity equation (∇·u = 0) is satisfied
    fn verify_continuity(&self, x: T, y: T, t: T) -> T {
        let vel = self.exact_velocity(x, y, t);
        let du_dx = self.velocity_derivative_x(x, y, t);
        let dv_dy = self.velocity_derivative_y(x, y, t);
        du_dx + dv_dy
    }

    /// Velocity derivatives (needed for source term computation)
    fn velocity_derivative_x(&self, x: T, y: T, t: T) -> T;
    fn velocity_derivative_y(&self, x: T, y: T, t: T) -> T;
    fn velocity_derivative_t(&self, x: T, y: T, t: T) -> Vector2<T>;

    /// Laplacian of velocity field
    fn velocity_laplacian(&self, x: T, y: T, t: T) -> Vector2<T>;
}

/// Manufactured solution for 2D incompressible Navier-Stokes using polynomial functions
/// This provides a complete MMS with analytical source terms
pub struct PolynomialNavierStokesMMS<T: RealField + Copy> {
    /// Kinematic viscosity
    pub nu: T,
    /// Density
    pub rho: T,
    /// Amplitude coefficients for velocity
    pub u_amp: T,
    pub v_amp: T,
    /// Amplitude coefficient for pressure
    pub p_amp: T,
}

impl<T: RealField + Copy + FromPrimitive> PolynomialNavierStokesMMS<T> {
    /// Create new polynomial MMS with specified parameters
    pub fn new(nu: T, rho: T, u_amp: T, v_amp: T, p_amp: T) -> Self {
        Self {
            nu,
            rho,
            u_amp,
            v_amp,
            p_amp,
        }
    }

    /// Create with default amplitudes for standard verification
    pub fn default(nu: T, rho: T) -> Self {
        Self::new(
            nu,
            rho,
            <T as FromPrimitive>::from_f64(1.0).unwrap(),
            <T as FromPrimitive>::from_f64(0.5).unwrap(),
            <T as FromPrimitive>::from_f64(0.1).unwrap(),
        )
    }
}

impl<T: RealField + Copy + FromPrimitive> NavierStokesManufacturedSolution<T>
    for PolynomialNavierStokesMMS<T>
{
    /// Exact velocity solution: u = A*sin(πx)*cos(πy)*exp(-2νπ²t)
    ///                       v = B*cos(πx)*sin(πy)*exp(-2νπ²t)
    fn exact_velocity(&self, x: T, y: T, t: T) -> Vector2<T> {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi * t,
        );

        let u = self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let v = self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        Vector2::new(u, v)
    }

    /// Exact pressure solution: p = C*sin(2πx)*cos(2πy)*exp(-4νπ²t)
    fn exact_pressure(&self, x: T, y: T, t: T) -> T {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(4.0).unwrap() * self.nu * pi * pi * t,
        );

        self.p_amp
            * nalgebra::ComplexField::sin(<T as FromPrimitive>::from_f64(2.0).unwrap() * pi * x)
            * nalgebra::ComplexField::cos(<T as FromPrimitive>::from_f64(2.0).unwrap() * pi * y)
            * decay
    }

    /// Source term for u-momentum equation: ∂u/∂t + u·∇u = -∇p/ρ + ν∇²u + f_u
    fn momentum_source_u(&self, x: T, y: T, t: T) -> T {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let two = <T as FromPrimitive>::from_f64(2.0).unwrap();
        let four = <T as FromPrimitive>::from_f64(4.0).unwrap();

        let decay = nalgebra::ComplexField::exp(-two * self.nu * pi * pi * t);
        let decay_4nu = nalgebra::ComplexField::exp(-four * self.nu * pi * pi * t);

        // ∂u/∂t
        let du_dt = -two
            * self.nu
            * pi
            * pi
            * self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;

        // u·∇u = u*∂u/∂x + v*∂u/∂y
        let u = self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let v = self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        let du_dx = self.u_amp
            * pi
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let du_dy = -self.u_amp
            * pi
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        let convection = u * du_dx + v * du_dy;

        // -∇p/ρ (pressure gradient contribution to u-momentum)
        let dp_dx = two
            * pi
            * self.p_amp
            * nalgebra::ComplexField::cos(two * pi * x)
            * nalgebra::ComplexField::cos(two * pi * y)
            * decay_4nu;
        let pressure_term = -dp_dx / self.rho;

        // ν∇²u
        let d2u_dx2 = -pi
            * pi
            * self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let d2u_dy2 = -pi
            * pi
            * self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let diffusion = self.nu * (d2u_dx2 + d2u_dy2);

        // Source term: f_u = ∂u/∂t + u·∇u + ∇p/ρ - ν∇²u
        du_dt + convection + pressure_term - diffusion
    }

    /// Source term for v-momentum equation
    fn momentum_source_v(&self, x: T, y: T, t: T) -> T {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let two = <T as FromPrimitive>::from_f64(2.0).unwrap();
        let four = <T as FromPrimitive>::from_f64(4.0).unwrap();

        let decay = nalgebra::ComplexField::exp(-two * self.nu * pi * pi * t);
        let decay_4nu = nalgebra::ComplexField::exp(-four * self.nu * pi * pi * t);

        // ∂v/∂t
        let dv_dt = -two
            * self.nu
            * pi
            * pi
            * self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        // u·∇v = u*∂v/∂x + v*∂v/∂y
        let u = self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;
        let v = self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        let dv_dx = -self.v_amp
            * pi
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;
        let dv_dy = self.v_amp
            * pi
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;

        let convection = u * dv_dx + v * dv_dy;

        // -∇p/ρ (pressure gradient contribution to v-momentum)
        let dp_dy = -two
            * pi
            * self.p_amp
            * nalgebra::ComplexField::sin(two * pi * x)
            * nalgebra::ComplexField::sin(two * pi * y)
            * decay_4nu;
        let pressure_term = -dp_dy / self.rho;

        // ν∇²v
        let d2v_dx2 = -pi
            * pi
            * self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;
        let d2v_dy2 = -pi
            * pi
            * self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;
        let diffusion = self.nu * (d2v_dx2 + d2v_dy2);

        // Source term: f_v = ∂v/∂t + u·∇v + ∇p/ρ - ν∇²v
        dv_dt + convection + pressure_term - diffusion
    }

    fn velocity_derivative_x(&self, x: T, y: T, t: T) -> T {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi * t,
        );
        self.u_amp
            * pi
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay
    }

    fn velocity_derivative_y(&self, x: T, y: T, t: T) -> T {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi * t,
        );
        -self.u_amp
            * pi
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay
    }

    fn velocity_derivative_t(&self, x: T, y: T, t: T) -> Vector2<T> {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let decay_factor = -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi;
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi * t,
        );

        let du_dt = self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay
            * decay_factor;
        let dv_dt = self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay
            * decay_factor;

        Vector2::new(du_dt, dv_dt)
    }

    fn velocity_laplacian(&self, x: T, y: T, t: T) -> Vector2<T> {
        let pi = <T as FromPrimitive>::from_f64(PI).unwrap();
        let pi_sq = pi * pi;
        let decay = nalgebra::ComplexField::exp(
            -<T as FromPrimitive>::from_f64(2.0).unwrap() * self.nu * pi * pi * t,
        );

        // ∇²u = ∂²u/∂x² + ∂²u/∂y² = -π²u_amp*sin(πx)cos(πy)*decay (twice)
        let lapl_u = -<T as FromPrimitive>::from_f64(2.0).unwrap()
            * pi_sq
            * self.u_amp
            * nalgebra::ComplexField::sin(pi * x)
            * nalgebra::ComplexField::cos(pi * y)
            * decay;

        // ∇²v = ∂²v/∂x² + ∂²v/∂y² = -π²v_amp*cos(πx)sin(πy)*decay (twice)
        let lapl_v = -<T as FromPrimitive>::from_f64(2.0).unwrap()
            * pi_sq
            * self.v_amp
            * nalgebra::ComplexField::cos(pi * x)
            * nalgebra::ComplexField::sin(pi * y)
            * decay;

        Vector2::new(lapl_u, lapl_v)
    }
}

/// Taylor-Green vortex solution for 2D Navier-Stokes
#[derive(Clone, Copy)]
pub struct TaylorGreenManufactured<T: RealField + Copy> {
    /// Kinematic viscosity
    pub nu: T,
    /// Wave number
    pub k: T,
}

impl<T: RealField + Copy + FromPrimitive> TaylorGreenManufactured<T> {
    /// Create a new Taylor-Green manufactured solution
    pub fn new(nu: T) -> Self {
        let pi = T::from_f64_or_one(PI);
        Self { nu, k: pi }
    }

    /// Get velocity components
    pub fn velocity(&self, x: T, y: T, t: T) -> Vector2<T> {
        let decay =
            nalgebra::ComplexField::exp(-T::from_f64_or_one(2.0) * self.nu * self.k * self.k * t);
        let u = nalgebra::ComplexField::sin(self.k * x)
            * nalgebra::ComplexField::cos(self.k * y)
            * decay;
        let v = -nalgebra::ComplexField::cos(self.k * x)
            * nalgebra::ComplexField::sin(self.k * y)
            * decay;
        Vector2::new(u, v)
    }

    /// Get pressure field
    pub fn pressure(&self, x: T, y: T, t: T) -> T {
        let decay =
            nalgebra::ComplexField::exp(-T::from_f64_or_one(4.0) * self.nu * self.k * self.k * t);
        let quarter = T::from_f64_or_one(0.25);
        -quarter
            * (nalgebra::ComplexField::cos(T::from_f64_or_one(2.0) * self.k * x)
                + nalgebra::ComplexField::cos(T::from_f64_or_one(2.0) * self.k * y))
            * decay
    }

    /// Get vorticity
    pub fn vorticity(&self, x: T, y: T, t: T) -> T {
        let decay =
            nalgebra::ComplexField::exp(-T::from_f64_or_one(2.0) * self.nu * self.k * self.k * t);
        -T::from_f64_or_one(2.0)
            * self.k
            * nalgebra::ComplexField::sin(self.k * x)
            * nalgebra::ComplexField::sin(self.k * y)
            * decay
    }

    /// Get kinetic energy
    pub fn kinetic_energy(&self, t: T) -> T {
        let decay =
            nalgebra::ComplexField::exp(-T::from_f64_or_one(4.0) * self.nu * self.k * self.k * t);
        T::from_f64_or_one(0.25) * decay
    }

    /// Get enstrophy (integral of vorticity squared)
    pub fn enstrophy(&self, t: T) -> T {
        let decay =
            nalgebra::ComplexField::exp(-T::from_f64_or_one(4.0) * self.nu * self.k * self.k * t);
        self.k * self.k * decay
    }
}

/// Kovasznay flow - exact solution for 2D steady Navier-Stokes
pub struct KovasznayFlow<T: RealField + Copy> {
    /// Reynolds number
    pub re: T,
    /// Lambda parameter
    pub lambda: T,
}

impl<T: RealField + Copy + FromPrimitive> KovasznayFlow<T> {
    /// Create a new Kovasznay flow solution
    pub fn new(re: T) -> Self {
        let half_re = re / T::from_f64_or_one(2.0);
        let pi_val = T::from_f64_or_one(PI);
        let discriminant = nalgebra::ComplexField::sqrt(
            half_re * half_re + T::from_f64_or_one(4.0) * pi_val * pi_val,
        );
        let lambda = half_re - discriminant;
        Self { re, lambda }
    }

    /// Get velocity components
    pub fn velocity(&self, x: T, y: T) -> Vector2<T> {
        let pi = T::from_f64_or_one(PI);
        let two_pi = T::from_f64_or_one(2.0) * pi;

        let exp_lambda_x = nalgebra::ComplexField::exp(self.lambda * x);
        let u = T::one() - exp_lambda_x * nalgebra::ComplexField::cos(two_pi * y);
        let v = self.lambda / two_pi * exp_lambda_x * nalgebra::ComplexField::sin(two_pi * y);

        Vector2::new(u, v)
    }

    /// Get pressure field
    pub fn pressure(&self, x: T) -> T {
        let half = T::from_f64_or_one(0.5);
        half * (T::one() - nalgebra::ComplexField::exp(T::from_f64_or_one(2.0) * self.lambda * x))
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
        let expected_ratio = (-4.0f64 * 0.01 * PI * PI).exp();
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
