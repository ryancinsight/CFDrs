//! Manufactured solutions for Burgers equation
//!
//! Burgers equation: ∂u/∂t + u(∂u/∂x) = ν(∂²u/∂x²)
//!
//! This equation combines nonlinear advection with diffusion, making it an excellent
//! test case for shock-capturing schemes and nonlinear convection-diffusion problems.
//!
//! References:
//! - Burgers, J.M. (1948) "A mathematical model illustrating the theory of turbulence"
//! - Bateman, H. (1915) "Some recent researches on the motion of fluids"
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"

use super::ManufacturedSolution;
use nalgebra::RealField;
use num_traits::Float;
use std::f64::consts::PI;

/// Manufactured solution for Burgers equation
///
/// Solution form: u(x,t) = a + b*sin(k*x - ω*t)
/// where parameters are chosen to satisfy Burgers equation with a source term
pub struct ManufacturedBurgers<T: RealField + Float> {
    /// Mean velocity
    a: T,
    /// Amplitude of sinusoidal component
    b: T,
    /// Wave number
    k: T,
    /// Angular frequency
    omega: T,
    /// Kinematic viscosity
    nu: T,
}

impl<T: RealField + Float> ManufacturedBurgers<T> {
    /// Create new manufactured Burgers solution
    ///
    /// # Arguments
    ///
    /// * `a` - Mean velocity (baseline value)
    /// * `b` - Amplitude of oscillation
    /// * `k` - Wave number (2π/wavelength)
    /// * `omega` - Angular frequency (2π/period)
    /// * `nu` - Kinematic viscosity
    ///
    /// # Example
    ///
    /// ```ignore
    /// use cfd_validation::manufactured::ManufacturedBurgers;
    ///
    /// // Create solution with viscosity 0.01
    /// let burgers = ManufacturedBurgers::new(1.0, 0.5, 2.0*std::f64::consts::PI, 1.0, 0.01);
    /// ```
    pub fn new(a: T, b: T, k: T, omega: T, nu: T) -> Self {
        Self { a, b, k, omega, nu }
    }

    /// Create default manufactured Burgers solution
    ///
    /// Uses standard parameters: a=1, b=0.5, k=2π, ω=1, ν=0.01
    pub fn default_solution() -> Self {
        let pi = T::from(PI).unwrap();
        let two = T::from(2.0).unwrap();
        Self::new(
            T::one(),                    // a = 1
            T::from(0.5).unwrap(),       // b = 0.5
            two * pi,                    // k = 2π
            T::one(),                    // ω = 1
            T::from(0.01).unwrap(),      // ν = 0.01
        )
    }
}

impl<T: RealField + Float> ManufacturedSolution<T> for ManufacturedBurgers<T> {
    /// Exact solution: u(x,t) = a + b*sin(k*x - ω*t)
    fn exact_solution(&self, x: T, _y: T, _z: T, t: T) -> T {
        let phase = self.k * x - self.omega * t;
        self.a + self.b * Float::sin(phase)
    }

    /// Source term required to satisfy Burgers equation
    ///
    /// For u = a + b*sin(θ) where θ = kx - ωt:
    ///
    /// ∂u/∂t = -bω*cos(θ)
    /// ∂u/∂x = bk*cos(θ)
    /// ∂²u/∂x² = -bk²*sin(θ)
    ///
    /// Burgers: ∂u/∂t + u(∂u/∂x) = ν(∂²u/∂x²)
    ///
    /// Source = ∂u/∂t + u(∂u/∂x) - ν(∂²u/∂x²)
    ///        = -bω*cos(θ) + (a + b*sin(θ))(bk*cos(θ)) + νbk²*sin(θ)
    ///        = -bω*cos(θ) + abk*cos(θ) + b²k*sin(θ)*cos(θ) + νbk²*sin(θ)
    fn source_term(&self, x: T, _y: T, _z: T, t: T) -> T {
        let phase = self.k * x - self.omega * t;
        let sin_theta = Float::sin(phase);
        let cos_theta = Float::cos(phase);

        // Temporal derivative term
        let temporal = -self.b * self.omega * cos_theta;

        // Nonlinear advection term: u(∂u/∂x)
        let u = self.a + self.b * sin_theta;
        let du_dx = self.b * self.k * cos_theta;
        let advection = u * du_dx;

        // Viscous diffusion term: -ν(∂²u/∂x²)
        let d2u_dx2 = -self.b * self.k * self.k * sin_theta;
        let viscous = -self.nu * d2u_dx2;

        temporal + advection + viscous
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_burgers_solution_properties() {
        let burgers = ManufacturedBurgers::<f64>::default_solution();

        // Test at origin and t=0
        let u_0 = burgers.exact_solution(0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(u_0, 1.0, epsilon = 1e-10); // a + b*sin(0) = 1 + 0 = 1

        // Test periodicity in space (period = 2π/k = 1 for k=2π)
        let x1 = 0.0;
        let x2 = 1.0; // One period
        let t = 0.5;
        let u1 = burgers.exact_solution(x1, 0.0, 0.0, t);
        let u2 = burgers.exact_solution(x2, 0.0, 0.0, t);
        assert_relative_eq!(u1, u2, epsilon = 1e-10);
    }

    #[test]
    fn test_burgers_source_term_nonzero() {
        let burgers = ManufacturedBurgers::<f64>::default_solution();

        // Source term should be non-zero for manufactured solution
        let source = burgers.source_term(0.5, 0.0, 0.0, 0.1);
        assert!(source.abs() > 1e-10, "Source term should be non-zero");
    }

    #[test]
    fn test_burgers_temporal_evolution() {
        let burgers = ManufacturedBurgers::<f64>::default_solution();

        let x = 0.5;
        let t1 = 0.0;
        let t2 = 0.1;

        let u1 = burgers.exact_solution(x, 0.0, 0.0, t1);
        let u2 = burgers.exact_solution(x, 0.0, 0.0, t2);

        // Solution should evolve in time
        assert!(
            (u1 - u2).abs() > 1e-10,
            "Solution should change over time"
        );
    }

    #[test]
    fn test_burgers_bounds() {
        let burgers = ManufacturedBurgers::<f64>::default_solution();

        // Solution should remain bounded: a-b <= u <= a+b
        // For a=1, b=0.5: 0.5 <= u <= 1.5
        for i in 0..100 {
            let x = f64::from(i) * 0.01;
            let t = 0.5;
            let u = burgers.exact_solution(x, 0.0, 0.0, t);
            assert!((0.5 - 1e-10..=1.5 + 1e-10).contains(&u), "u={u} out of bounds");
        }
    }
}
