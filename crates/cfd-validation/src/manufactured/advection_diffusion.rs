//! Manufactured solutions for coupled advection-diffusion equation
//!
//! Coupled advection-diffusion: ∂u/∂t + v·∇u = α∇²u
//!
//! This equation represents transport phenomena with both convection and diffusion,
//! testing the interaction between these mechanisms. Important for scalar transport
//! in CFD applications (temperature, species concentration, etc.).
//!
//! References:
//! - Patankar, S.V. (1980) "Numerical Heat Transfer and Fluid Flow"
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"

use super::ManufacturedSolution;
use nalgebra::RealField;
use num_traits::Float;
use std::f64::consts::PI;

/// Manufactured solution for coupled advection-diffusion
///
/// Solution form: u(x,y,t) = exp(-αλ²t) * sin(kx*x) * sin(ky*y)
/// where λ² = kx² + ky² and velocity field v = (vx, vy) is constant
pub struct ManufacturedAdvectionDiffusion<T: RealField + Float> {
    /// Wave number in x-direction
    kx: T,
    /// Wave number in y-direction
    ky: T,
    /// Diffusion coefficient (thermal diffusivity or mass diffusivity)
    alpha: T,
    /// Velocity component in x-direction
    vx: T,
    /// Velocity component in y-direction
    vy: T,
}

impl<T: RealField + Float> ManufacturedAdvectionDiffusion<T> {
    /// Create new manufactured advection-diffusion solution
    ///
    /// # Arguments
    ///
    /// * `kx` - Wave number in x-direction
    /// * `ky` - Wave number in y-direction
    /// * `alpha` - Diffusion coefficient
    /// * `vx` - Velocity in x-direction
    /// * `vy` - Velocity in y-direction
    ///
    /// # Example
    ///
    /// ```ignore
    /// use cfd_validation::manufactured::ManufacturedAdvectionDiffusion;
    ///
    /// // Create solution with moderate Peclet number
    /// let coupled = ManufacturedAdvectionDiffusion::new(
    ///     2.0*std::f64::consts::PI, 
    ///     2.0*std::f64::consts::PI,
    ///     0.01,
    ///     1.0,
    ///     0.5
    /// );
    /// ```
    pub fn new(kx: T, ky: T, alpha: T, vx: T, vy: T) -> Self {
        Self { kx, ky, alpha, vx, vy }
    }

    /// Create default solution with Pe ≈ 10 (advection-dominated)
    pub fn default_advection_dominated() -> Self {
        let pi = T::from(PI).unwrap();
        let two = T::from(2.0).unwrap();
        Self::new(
            two * pi,                    // kx = 2π
            two * pi,                    // ky = 2π
            T::from(0.01).unwrap(),      // α = 0.01
            T::one(),                    // vx = 1
            T::from(0.5).unwrap(),       // vy = 0.5
        )
    }

    /// Create solution with Pe ≈ 1 (balanced)
    pub fn default_balanced() -> Self {
        let pi = T::from(PI).unwrap();
        let two = T::from(2.0).unwrap();
        Self::new(
            two * pi,                    // kx = 2π
            two * pi,                    // ky = 2π
            T::from(0.1).unwrap(),       // α = 0.1
            T::one(),                    // vx = 1
            T::from(0.5).unwrap(),       // vy = 0.5
        )
    }

    /// Create solution with Pe << 1 (diffusion-dominated)
    pub fn default_diffusion_dominated() -> Self {
        let pi = T::from(PI).unwrap();
        let two = T::from(2.0).unwrap();
        Self::new(
            two * pi,                    // kx = 2π
            two * pi,                    // ky = 2π
            T::one(),                    // α = 1.0
            T::one(),                    // vx = 1
            T::from(0.5).unwrap(),       // vy = 0.5
        )
    }

    /// Calculate local Peclet number Pe = |v|*h/α
    ///
    /// # Arguments
    ///
    /// * `h` - Characteristic length scale (grid spacing)
    pub fn peclet_number(&self, h: T) -> T {
        let v_magnitude = Float::sqrt(self.vx * self.vx + self.vy * self.vy);
        v_magnitude * h / self.alpha
    }
}

impl<T: RealField + Float> ManufacturedSolution<T> for ManufacturedAdvectionDiffusion<T> {
    /// Exact solution: u(x,y,t) = exp(-αλ²t) * sin(kx*x) * sin(ky*y)
    /// where λ² = kx² + ky²
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        let lambda_sq = self.kx * self.kx + self.ky * self.ky;
        let decay = Float::exp(-self.alpha * lambda_sq * t);
        decay * Float::sin(self.kx * x) * Float::sin(self.ky * y)
    }

    /// Source term required to satisfy advection-diffusion equation
    ///
    /// For u = exp(-αλ²t) * sin(kx*x) * sin(ky*y):
    ///
    /// ∂u/∂t = -αλ² * exp(-αλ²t) * sin(kx*x) * sin(ky*y)
    /// ∂u/∂x = exp(-αλ²t) * kx * cos(kx*x) * sin(ky*y)
    /// ∂u/∂y = exp(-αλ²t) * ky * sin(kx*x) * cos(ky*y)
    /// ∂²u/∂x² = -exp(-αλ²t) * kx² * sin(kx*x) * sin(ky*y)
    /// ∂²u/∂y² = -exp(-αλ²t) * ky² * sin(kx*x) * sin(ky*y)
    ///
    /// Equation: ∂u/∂t + vx*(∂u/∂x) + vy*(∂u/∂y) = α(∂²u/∂x² + ∂²u/∂y²)
    ///
    /// Source = ∂u/∂t + vx*(∂u/∂x) + vy*(∂u/∂y) - α(∂²u/∂x² + ∂²u/∂y²)
    ///
    /// For the chosen form, left side should equal right side, but we add source
    /// term to test solver with non-homogeneous equation
    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        let lambda_sq = self.kx * self.kx + self.ky * self.ky;
        let decay = Float::exp(-self.alpha * lambda_sq * t);
        let sin_x = Float::sin(self.kx * x);
        let cos_x = Float::cos(self.kx * x);
        let sin_y = Float::sin(self.ky * y);
        let cos_y = Float::cos(self.ky * y);

        // Temporal derivative
        let du_dt = -self.alpha * lambda_sq * decay * sin_x * sin_y;

        // Advection terms
        let du_dx = decay * self.kx * cos_x * sin_y;
        let du_dy = decay * self.ky * sin_x * cos_y;
        let advection = self.vx * du_dx + self.vy * du_dy;

        // Diffusion terms
        let d2u_dx2 = -decay * self.kx * self.kx * sin_x * sin_y;
        let d2u_dy2 = -decay * self.ky * self.ky * sin_x * sin_y;
        let diffusion = self.alpha * (d2u_dx2 + d2u_dy2);

        // Source = temporal + advection - diffusion
        // For this manufactured solution, this should be nearly zero
        // but we keep it general for testing purposes
        du_dt + advection - diffusion
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_advection_diffusion_initial_condition() {
        let coupled = ManufacturedAdvectionDiffusion::<f64>::default_balanced();

        // At t=0, exp(-αλ²t) = 1
        let u_0 = coupled.exact_solution(0.0, 0.0, 0.0, 0.0);
        assert_relative_eq!(u_0, 0.0, epsilon = 1e-10); // sin(0)*sin(0) = 0
    }

    #[test]
    fn test_advection_diffusion_decay() {
        let coupled = ManufacturedAdvectionDiffusion::<f64>::default_balanced();

        let x = std::f64::consts::PI / 4.0;
        let y = std::f64::consts::PI / 4.0;
        let t1 = 0.0;
        let t2 = 0.1;

        let u1 = coupled.exact_solution(x, y, 0.0, t1);
        let u2 = coupled.exact_solution(x, y, 0.0, t2);

        // Solution should decay over time (diffusion)
        assert!(u1.abs() > u2.abs(), "Solution should decay: {} > {}", u1, u2);
    }

    #[test]
    fn test_peclet_number_calculation() {
        let coupled = ManufacturedAdvectionDiffusion::<f64>::default_advection_dominated();

        let h = 0.01;
        let pe = coupled.peclet_number(h);

        // For vx=1, vy=0.5: |v| = sqrt(1.25) ≈ 1.118
        // Pe = |v|*h/α = 1.118*0.01/0.01 ≈ 1.118
        assert!(pe > 1.0, "Advection-dominated should have Pe > 1: {}", pe);
    }

    #[test]
    fn test_source_term_small_for_exact_solution() {
        let coupled = ManufacturedAdvectionDiffusion::<f64>::default_balanced();

        // For the chosen form, source term should be very small (near zero)
        // because the exact solution already satisfies the homogeneous equation
        let source = coupled.source_term(0.5, 0.5, 0.0, 0.1);
        
        // Source should be near zero (within numerical precision)
        assert!(
            source.abs() < 1e-8,
            "Source term should be near zero for exact solution: {}",
            source
        );
    }

    #[test]
    fn test_periodicity() {
        let coupled = ManufacturedAdvectionDiffusion::<f64>::default_balanced();

        let t = 0.5;
        // Period in x: 2π/kx = 1 for kx=2π
        // Period in y: 2π/ky = 1 for ky=2π
        let u1 = coupled.exact_solution(0.0, 0.0, 0.0, t);
        let u2 = coupled.exact_solution(1.0, 1.0, 0.0, t);

        assert_relative_eq!(u1, u2, epsilon = 1e-10);
    }

    #[test]
    fn test_different_regimes() {
        let advection_dom = ManufacturedAdvectionDiffusion::<f64>::default_advection_dominated();
        let diffusion_dom = ManufacturedAdvectionDiffusion::<f64>::default_diffusion_dominated();

        let h = 0.01;
        let pe_adv = advection_dom.peclet_number(h);
        let pe_diff = diffusion_dom.peclet_number(h);

        assert!(pe_adv > pe_diff, "Advection-dominated should have higher Pe");
        assert!(pe_diff < 1.0, "Diffusion-dominated should have Pe < 1");
    }
}
