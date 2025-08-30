//! Courant-Friedrichs-Lewy (CFL) stability conditions for numerical schemes
//!
//! References:
//! - Courant, Friedrichs, Lewy (1928) "Über die partiellen Differenzengleichungen der mathematischen Physik"
//! - Hirsch (1988) "Numerical Computation of Internal and External Flows"

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// CFL condition calculator for various numerical schemes
pub struct CFLCalculator<T: RealField + Copy> {
    /// Grid spacing in x-direction
    dx: T,
    /// Grid spacing in y-direction
    dy: T,
    /// Time step
    dt: T,
}

impl<T: RealField + Copy + FromPrimitive> CFLCalculator<T> {
    /// Create new CFL calculator
    pub fn new(dx: T, dy: T, dt: T) -> Self {
        Self { dx, dy, dt }
    }

    /// Calculate CFL number for pure advection
    /// CFL = u*dt/dx + v*dt/dy
    ///
    /// Stability limits:
    /// - Explicit Euler: CFL ≤ 1.0
    /// - Lax-Wendroff: CFL ≤ 1.0
    /// - MacCormack: CFL ≤ 1.0
    /// - QUICK: CFL ≤ 0.75 (for third-order accuracy)
    pub fn advection_cfl(&self, u: T, v: T) -> T {
        u.abs() * self.dt / self.dx + v.abs() * self.dt / self.dy
    }

    /// Calculate diffusion number (von Neumann number) for pure diffusion
    /// D = ν*dt/dx² + ν*dt/dy²
    ///
    /// Stability limits:
    /// - Explicit Euler: D ≤ 0.5 (2D), D ≤ 0.25 (3D)
    /// - Crank-Nicolson: Unconditionally stable
    pub fn diffusion_number(&self, nu: T) -> T {
        nu * self.dt * (T::one() / (self.dx * self.dx) + T::one() / (self.dy * self.dy))
    }

    /// Calculate combined CFL for advection-diffusion problems
    ///
    /// For explicit schemes, both conditions must be satisfied:
    /// - CFL ≤ CFL_max (advection stability)
    /// - D ≤ D_max (diffusion stability)
    ///
    /// The Peclet number Pe = u*dx/ν determines relative importance
    pub fn combined_stability(&self, u: T, v: T, nu: T) -> (T, T, T) {
        let cfl = self.advection_cfl(u, v);
        let diff = self.diffusion_number(nu);
        let peclet = u.abs() * self.dx / nu;
        (cfl, diff, peclet)
    }

    /// Check stability for explicit Euler time stepping
    pub fn is_stable_explicit_euler(&self, u: T, v: T, nu: T) -> bool {
        let cfl = self.advection_cfl(u, v);
        let diff = self.diffusion_number(nu);

        // Explicit Euler requires CFL ≤ 1 and D ≤ 0.5 in 2D
        cfl <= T::one() && diff <= T::from_f64(0.5).unwrap_or_else(T::one)
    }

    /// Check stability for QUICK scheme with explicit time stepping
    pub fn is_stable_quick(&self, u: T, v: T, nu: T) -> bool {
        let cfl = self.advection_cfl(u, v);
        let diff = self.diffusion_number(nu);

        // QUICK requires stricter CFL for third-order accuracy
        // Leonard (1979) suggests CFL ≤ 0.75 for stability
        let cfl_limit = T::from_f64(0.75).unwrap_or_else(T::one);
        let diff_limit = T::from_f64(0.5).unwrap_or_else(T::one);

        cfl <= cfl_limit && diff <= diff_limit
    }

    /// Calculate maximum stable time step for given flow conditions
    pub fn max_stable_dt(&self, u_max: T, v_max: T, nu: T) -> T {
        // CFL constraint: dt ≤ min(dx/|u|, dy/|v|)
        let dt_advection = if u_max > T::zero() && v_max > T::zero() {
            let dt_x = self.dx / u_max.abs();
            let dt_y = self.dy / v_max.abs();
            dt_x.min(dt_y)
        } else if u_max > T::zero() {
            self.dx / u_max.abs()
        } else if v_max > T::zero() {
            self.dy / v_max.abs()
        } else {
            T::from_f64(1e10).unwrap_or_else(T::one) // Large value for zero velocity
        };

        // Diffusion constraint: dt ≤ 0.5 * min(dx², dy²) / ν
        let dt_diffusion = if nu > T::zero() {
            let dx2 = self.dx * self.dx;
            let dy2 = self.dy * self.dy;
            T::from_f64(0.5).unwrap_or_else(T::one) * dx2.min(dy2) / nu
        } else {
            T::from_f64(1e10).unwrap_or_else(T::one) // Large value for zero diffusion
        };

        // Return minimum of both constraints
        dt_advection.min(dt_diffusion)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cfl_calculation() {
        let calculator = CFLCalculator::new(0.01, 0.01, 0.001);

        // Test pure advection
        let cfl = calculator.advection_cfl(1.0, 0.5);
        assert!((cfl - 0.15).abs() < 1e-10);

        // Test diffusion number
        let diff = calculator.diffusion_number(0.01);
        assert!((diff - 0.02).abs() < 1e-10);
    }

    #[test]
    fn test_stability_checks() {
        let calculator = CFLCalculator::new(0.01, 0.01, 0.001);

        // Should be stable for small velocities
        assert!(calculator.is_stable_explicit_euler(0.5, 0.5, 0.01));

        // Should be unstable for large velocities
        assert!(!calculator.is_stable_explicit_euler(20.0, 20.0, 0.01));
    }
}
