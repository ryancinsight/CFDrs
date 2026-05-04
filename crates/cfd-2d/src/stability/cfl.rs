//! Courant-Friedrichs-Lewy (CFL) stability conditions for numerical schemes
//!
//! References:
//! - Courant, Friedrichs, Lewy (1928) "Über die partiellen Differenzengleichungen der mathematischen Physik"
//! - Hirsch (1988) "Numerical Computation of Internal and External Flows"
//!
//! # Theorem (CFL Necessary Condition — Courant, Friedrichs, Lewy 1928)
//!
//! For an explicit time-marching scheme applied to the advection equation
//! $\partial u/\partial t + c\,\partial u/\partial x = 0$, stability requires
//! $\text{CFL} = |c|\Delta t / \Delta x \le C_{\max}$, where $C_{\max} \le 1$
//! depends on the spatial discretisation.
//!
//! **Proof sketch**:
//! The domain of dependence of the PDE at $(x, t+\Delta t)$ is the characteristic
//! interval $[x - c\Delta t,\, x + c\Delta t]$. The numerical domain of dependence
//! spans $[x - \Delta x,\, x + \Delta x]$ for a 3-point stencil. If the physical
//! domain of dependence is not contained in the numerical one, the scheme cannot
//! represent the correct solution. For 2D with velocities $(u, v)$ and diffusivity
//! $\nu$, the combined condition is $|u|\Delta t/\Delta x + |v|\Delta t/\Delta y
//! + 2\nu\Delta t(1/\Delta x^2 + 1/\Delta y^2) \le 1$.

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
    /// - `MacCormack`: CFL ≤ 1.0
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
    /// - CFL ≤ `CFL_max` (advection stability)
    /// - D ≤ `D_max` (diffusion stability)
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
        cfl <= T::one() && diff <= T::from_f64(0.5).expect("Exact mathematically representable f64")
    }

    /// Check stability for QUICK scheme with explicit time stepping
    pub fn is_stable_quick(&self, u: T, v: T, nu: T) -> bool {
        let cfl = self.advection_cfl(u, v);
        let diff = self.diffusion_number(nu);

        // QUICK requires stricter CFL for third-order accuracy
        // Leonard (1979) suggests CFL ≤ 0.75 for stability
        let cfl_limit = T::from_f64(0.75).expect("Exact mathematically representable f64");
        let diff_limit = T::from_f64(0.5).expect("Exact mathematically representable f64");

        cfl <= cfl_limit && diff <= diff_limit
    }

    /// Calculate maximum stable time step for given flow conditions.
    ///
    /// # Theorem
    /// For an unsplit explicit 2D advection-diffusion update, stability requires
    /// the summed advective and diffusive rates to remain bounded:
    ///
    /// ```text
    /// Δt_adv ≤ 1 / (|u|/Δx + |v|/Δy)
    /// Δt_diff ≤ 0.5 / (ν(1/Δx² + 1/Δy²))
    /// ```
    ///
    /// **Proof sketch**: the documented CFL condition is
    /// `|u|Δt/Δx + |v|Δt/Δy ≤ 1`; solving for `Δt` gives the reciprocal
    /// summed advective rate. The explicit 2D diffusion von Neumann condition is
    /// `νΔt(1/Δx² + 1/Δy²) ≤ 1/2`; solving for `Δt` gives the reciprocal
    /// summed diffusive rate. Taking the minimum satisfies both constraints.
    pub fn max_stable_dt(&self, u_max: T, v_max: T, nu: T) -> T {
        let advective_rate = u_max.abs() / self.dx + v_max.abs() / self.dy;
        let large_dt = T::from_f64(1e10).expect("Exact mathematically representable f64");

        let dt_advection = if advective_rate > T::zero() {
            T::one() / advective_rate
        } else {
            large_dt
        };

        // Diffusion constraint: dt ≤ 0.5 / (ν(1/dx² + 1/dy²)).
        let dt_diffusion = if nu > T::zero() {
            let inv_dx2 = T::one() / (self.dx * self.dx);
            let inv_dy2 = T::one() / (self.dy * self.dy);
            T::from_f64(0.5).expect("Exact mathematically representable f64")
                / (nu * (inv_dx2 + inv_dy2))
        } else {
            large_dt
        };

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
        assert!((cfl - 0.15_f64).abs() < 1e-10);

        // Test diffusion number
        let diff = calculator.diffusion_number(0.01);
        assert!((diff - 0.2_f64).abs() < 1e-10);
    }

    #[test]
    fn test_stability_checks() {
        let calculator = CFLCalculator::new(0.01, 0.01, 0.001);

        // Should be stable for small velocities
        assert!(calculator.is_stable_explicit_euler(0.5, 0.5, 0.01));

        // Should be unstable for large velocities
        assert!(!calculator.is_stable_explicit_euler(20.0, 20.0, 0.01));
    }

    #[test]
    fn max_stable_dt_respects_summed_2d_advection_cfl() {
        let calculator = CFLCalculator::new(0.01, 0.02, 0.001);
        let dt = calculator.max_stable_dt(2.0, 4.0, 0.0);
        let expected: f64 = 1.0 / (2.0 / 0.01 + 4.0 / 0.02);

        assert!((dt - expected).abs() < 1e-15);

        let stable_calculator = CFLCalculator::new(0.01, 0.02, dt);
        assert!(
            stable_calculator.advection_cfl(2.0, 4.0) <= 1.0 + 1e-14,
            "max_stable_dt must satisfy the summed 2D advection CFL"
        );
    }

    #[test]
    fn max_stable_dt_respects_summed_2d_diffusion_bound() {
        let calculator = CFLCalculator::new(0.01, 0.02, 0.001);
        let nu = 1.5e-5;
        let dt = calculator.max_stable_dt(0.0, 0.0, nu);
        let expected: f64 = 0.5 / (nu * (1.0 / 0.01_f64.powi(2) + 1.0 / 0.02_f64.powi(2)));

        assert!((dt - expected).abs() < 1e-15);

        let stable_calculator = CFLCalculator::new(0.01, 0.02, dt);
        assert!(
            stable_calculator.diffusion_number(nu) <= 0.5 + 1e-14,
            "max_stable_dt must satisfy the 2D diffusion von Neumann bound"
        );
    }
}
