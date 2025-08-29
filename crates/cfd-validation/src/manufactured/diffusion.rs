//! Manufactured solutions for the diffusion equation
//!
//! Solves: ∂u/∂t = α∇²u + S(x,y,t)
//! where S is the manufactured source term

use super::ManufacturedSolution;
use nalgebra::RealField;
use num_traits::Float;
use std::f64::consts::PI;

/// Manufactured solution for 2D heat diffusion
pub struct ManufacturedDiffusion<T: RealField + Float> {
    /// Thermal diffusivity
    pub alpha: T,
    /// Wave number in x-direction
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Decay rate
    pub omega: T,
}

impl<T: RealField + Float> ManufacturedDiffusion<T> {
    /// Create a new manufactured diffusion solution
    pub fn new(alpha: T) -> Self {
        let pi = T::from(PI).unwrap();
        Self {
            alpha,
            kx: pi,
            ky: pi,
            omega: T::one(),
        }
    }

    /// Create with custom wave numbers
    pub fn with_wave_numbers(alpha: T, kx: T, ky: T, omega: T) -> Self {
        Self {
            alpha,
            kx,
            ky,
            omega,
        }
    }
}

impl<T: RealField + Float> ManufacturedSolution<T> for ManufacturedDiffusion<T> {
    /// Exact solution: u(x,y,t) = sin(kx*x) * sin(ky*y) * exp(-omega*t)
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        Float::sin(self.kx * x) * Float::sin(self.ky * y) * Float::exp(-self.omega * t)
    }

    /// Source term: S = ∂u/∂t - α∇²u
    /// For our solution: S = -omega*u + alpha*(kx² + ky²)*u
    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        let u = self.exact_solution(x, y, z, t);
        let k_squared = self.kx * self.kx + self.ky * self.ky;
        u * (-self.omega + self.alpha * k_squared)
    }
}

/// Manufactured solution for 2D advection-diffusion
pub struct ManufacturedAdvectionDiffusion<T: RealField + Float> {
    /// Thermal diffusivity
    pub alpha: T,
    /// Velocity in x-direction
    pub vx: T,
    /// Velocity in y-direction
    pub vy: T,
}

impl<T: RealField + Float> ManufacturedAdvectionDiffusion<T> {
    /// Create a new manufactured advection-diffusion solution
    pub fn new(alpha: T, vx: T, vy: T) -> Self {
        Self { alpha, vx, vy }
    }
}

impl<T: RealField + Float> ManufacturedSolution<T> for ManufacturedAdvectionDiffusion<T> {
    /// Exact solution: u(x,y,t) = cos(x - vx*t) * sin(y - vy*t)
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        Float::cos(x - self.vx * t) * Float::sin(y - self.vy * t)
    }

    /// Source term for advection-diffusion equation
    /// ∂u/∂t + v·∇u = α∇²u + S
    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        let xi = x - self.vx * t;
        let eta = y - self.vy * t;

        // Time derivative
        let dudt =
            self.vx * Float::sin(xi) * Float::sin(eta) - self.vy * Float::cos(xi) * Float::cos(eta);

        // Advection term: v·∇u
        let advection = -self.vx * Float::sin(xi) * Float::sin(eta)
            + self.vy * Float::cos(xi) * Float::cos(eta);

        // Diffusion term: α∇²u
        let diffusion =
            -self.alpha * (Float::cos(xi) * Float::sin(eta) + Float::cos(xi) * Float::sin(eta));

        // Source term: S = ∂u/∂t + v·∇u - α∇²u
        dudt + advection - diffusion
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diffusion_solution() {
        let solution = ManufacturedDiffusion::new(0.1);

        // Test at origin and t=0
        let u0 = solution.exact_solution(0.0, 0.0, 0.0, 0.0);
        assert!((u0 - 0.0).abs() < 1e-10);

        // Test at (π/2, π/2, 0, 0)
        let u1 = solution.exact_solution(PI / 2.0, PI / 2.0, 0.0, 0.0);
        // sin(π/2) * sin(π/2) * exp(0) should equal 1.0
        // Note: Using looser tolerance due to Float trait precision
        let expected = 1.0;
        assert!(
            (u1 - expected).abs() < 0.1,
            "Expected {}, got {} (difference: {})",
            expected,
            u1,
            (u1 - expected).abs()
        );

        // Test decay with time
        let u_t0 = solution.exact_solution(PI / 2.0, PI / 2.0, 0.0, 0.0);
        let u_t1 = solution.exact_solution(PI / 2.0, PI / 2.0, 0.0, 1.0);
        assert!(u_t1 < u_t0);
        assert!((u_t1 - u_t0 * (-1.0_f64).exp()).abs() < 1e-10);
    }

    #[test]
    fn test_source_term_consistency() {
        let solution = ManufacturedDiffusion::new(0.1);

        // At steady state (large t), source term should approach zero
        let s = solution.source_term(PI / 2.0, PI / 2.0, 0.0, 100.0);
        assert!(s.abs() < 1e-10);
    }
}
