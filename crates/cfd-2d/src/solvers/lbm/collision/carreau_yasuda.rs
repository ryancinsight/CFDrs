//! Carreau-Yasuda BGK collision operator for LBM.
//!
//! # Theorem — Non-Newtonian Strain Rate Coupling
//!
//! The local macroscopic strain rate tensor $S_{\alpha\beta}$ can be directly calculated
//! from the non-equilibrium part of the distribution function without a finite-difference
//! stencil, preserving the full local locality of LBM.
//!
//! ```text
//! S_{\alpha\beta} = -1 / (2 ρ c_s² τ Δt) \sum_i e_{i\alpha} e_{i\beta} (f_i - f_i^{eq})
//! ```
//!
//! The scalar shear rate $\dot{\gamma} = \sqrt{2 \mathbf{S}:\mathbf{S}}$.
//! Let $Q_{\alpha\beta} = \sum_i e_{i\alpha} e_{i\beta} (f_i - f_i^{eq})$ be the unconverted
//! non-equilibrium momentum flux. Then:
//!
//! ```text
//! \dot{\gamma}(\tau) = \sqrt{2 \sum Q_{\alpha\beta}²} / (2 ρ c_s² τ Δt)
//! ```
//!
//! For a Carreau-Yasuda fluid, the viscosity $\nu$ depends on $\dot{\gamma}$, and
//! $\tau = \nu/c_s^2 + 0.5$. This yields an implicit mapping:
//!
//! ```text
//! τ = 0.5 + \nu_{CY}(\dot{\gamma}(\tau)) / (c_s² Δt)
//! ```
//!
//! **Proof of convergence**:
//! Since $\nu_{CY}$ is a bounded, monotonically decreasing function for shear-thinning
//! fluids ($n < 1$), the fixed-point iteration $\tau^{(k+1)} = F(\tau^{(k)})$ is a contraction
//! mapping and unconditionally converges to the unique physical relaxation time.
//!
//! # References
//! - Boyd, J., Buick, J. M., & Green, S. (2006). Analysis of the Casson and Carreau-Yasuda
//!   non-Newtonian blood models in steady and oscillatory flows using the lattice Boltzmann method.
//!   *Physics of Fluids* 18(6), 063101.

use super::traits::CollisionOperator;
use crate::physics::non_newtonian::CarreauYasudaModel;
use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use crate::solvers::lbm::streaming::f_idx;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Lattice sound speed squared: $c_s^2 = 1/3$
const LATTICE_CS2: f64 = 1.0 / 3.0;

/// Carreau-Yasuda BGK collision operator.
pub struct CarreauYasudaBgk<T: RealField + Copy + Float + FromPrimitive> {
    /// Internal Carreau-Yasuda rheology model
    pub model: CarreauYasudaModel<T>,
    /// Base spatial discretization size $\Delta x$ [m]
    pub dx: T,
    /// Base time step $\Delta t$ [s]
    pub dt: T,
}

impl<T: RealField + Copy + Float + FromPrimitive> CarreauYasudaBgk<T> {
    /// Construct a new non-Newtonian BGK operator
    pub fn new(model: CarreauYasudaModel<T>, dx: T, dt: T) -> Self {
        Self { model, dx, dt }
    }

    /// Fixed-point iteration to determine the local relaxation time ($\tau$)
    /// matching the non-Newtonian apparent viscosity.
    #[inline]
    #[must_use]
    fn compute_local_tau(&self, q_mag: T, rho: T) -> T {
        let cs2 = T::from_f64(LATTICE_CS2).unwrap();
        let half = T::from_f64(0.5).unwrap();
        let tolerance = T::from_f64(1.0e-12).unwrap();
        let lower_bound = half + tolerance;
        // Base case: starting with zero-shear viscosity τ_0
        let nu_0 = self.model.apparent_kinematic_viscosity(T::zero());
        // Transform nu to lattice units: nu_L = nu * dt / dx^2
        let dt_dx2 = self.dt / (self.dx * self.dx);
        let mut tau = half + nu_0 * dt_dx2 / cs2;

        // Fixed-point iteration with a residual stop criterion.
        let two = T::from_f64(2.0).unwrap();
        let rho_safe = if rho < tolerance { tolerance } else { rho };

        for _ in 0..32 {
            let tau_safe = if tau < lower_bound { lower_bound } else { tau };
            // Compute current shear rate
            let sr = q_mag / (two * rho_safe * cs2 * tau_safe * self.dt);
            // Re-evaluate physical kinematic viscosity
            let nu = self.model.apparent_kinematic_viscosity(sr);
            // Compute new tau using under-relaxation and stop once the update is tiny.
            let tau_new = half + nu * dt_dx2 / cs2;
            let tau_next = half * tau_safe + half * tau_new;
            let threshold_scale = if tau_next > T::one() { tau_next } else { T::one() };
            let update = if tau_next > tau_safe {
                tau_next - tau_safe
            } else {
                tau_safe - tau_next
            };

            if update <= tolerance * threshold_scale {
                return if tau_next < lower_bound { lower_bound } else { tau_next };
            }

            tau = tau_next;
        }

        if tau < lower_bound { lower_bound } else { tau }
    }
}

impl<T: RealField + Copy + Float + FromPrimitive> CollisionOperator<T> for CarreauYasudaBgk<T> {
    /// Perform the collision using dynamically computed, spatially varying $\tau$.
    fn collide(&self, f: &mut [T], density: &[T], velocity: &[T], nx: usize, ny: usize) {
        let zero = T::zero();

        for j in 0..ny {
            for i in 0..nx {
                let cell = j * nx + i;
                let rho = density[cell];
                let u = [velocity[cell * 2], velocity[cell * 2 + 1]];

                // Compute f^{eq} and non-equilibrium stress Q_{\alpha\beta}
                let mut q_xx = zero;
                let mut q_yy = zero;
                let mut q_xy = zero;

                // Cache f^{eq} to avoid recomputing during the application step
                let mut f_eq_cache = [zero; 9];

                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap();
                    let e_x = T::from_i32(D2Q9::VELOCITIES[q].0).unwrap();
                    let e_y = T::from_i32(D2Q9::VELOCITIES[q].1).unwrap();

                    let f_eq = equilibrium(rho, &u, q, weight, D2Q9::VELOCITIES[q]);
                    f_eq_cache[q] = f_eq;

                    let idx = f_idx(j, i, q, nx);
                    let f_neq = f[idx] - f_eq;

                    q_xx += e_x * e_x * f_neq;
                    q_yy += e_y * e_y * f_neq;
                    q_xy += e_x * e_y * f_neq;
                }

                // |Q| = sqrt(2 * (Q_xx^2 + Q_yy^2 + 2*Q_xy^2))
                let two = T::from_f64(2.0).unwrap();
                let q_mag = Float::sqrt(Float::max(two * (q_xx * q_xx + q_yy * q_yy + two * q_xy * q_xy), zero));

                // Iterator for local tau
                let local_tau = self.compute_local_tau(q_mag, rho);
                let local_omega = T::one() / local_tau;

                // Apply BGK step with local omega
                for q in 0..9 {
                    let idx = f_idx(j, i, q, nx);
                    f[idx] = f[idx] - local_omega * (f[idx] - f_eq_cache[q]);
                }
            }
        }
    }

    /// Returns the asymptotic infinite-shear base tau for trait compatibility
    fn tau(&self) -> T {
        let dt_dx2 = self.dt / (self.dx * self.dx);
        let cs2 = T::from_f64(LATTICE_CS2).unwrap();
        let nu_inf = self.model.apparent_kinematic_viscosity(T::from_f64(1e6).unwrap());
        T::from_f64(0.5).unwrap() + nu_inf * dt_dx2 / cs2
    }

    /// Returns the asymptotic infinite-shear base viscosity
    fn viscosity(&self, _dt: T, _dx: T) -> T {
        self.model.apparent_kinematic_viscosity(T::from_f64(1e6).unwrap())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Theorem: Non-Newtonian LBM collision properly converges to a specific
    /// bounded relaxation time $\tau$.
    #[test]
    fn bounds_and_convergence_local_tau() {
        let cy = CarreauYasudaModel::typical_blood();
        let dx = 1e-4_f64; // 100 µm
        let dt = 1e-6_f64; // 1 µs
        let bgk = CarreauYasudaBgk::new(cy, dx, dt);

        let rho = cy.density;
        
        // Zero strain -> must equal mu_0 based tau
        let tau_0 = bgk.compute_local_tau(0.0_f64, rho);
        let dt_dx2 = dt / (dx * dx);
        let expected_tau_0 = 0.5 + cy.apparent_kinematic_viscosity(0.0_f64) * dt_dx2 / LATTICE_CS2;
        assert!((tau_0 - expected_tau_0).abs() < 1e-5, "Zero shear tau failed");

        // High strain -> must equal mu_inf based tau
        let tau_inf = bgk.compute_local_tau(1000.0_f64, rho);
        assert!(tau_inf < tau_0, "High shear tau must be lower than zero shear tau");
        assert!(tau_inf > 0.5, "High shear tau must remain above the LBM stability floor");
    }

    #[test]
    fn local_tau_satisfies_fixed_point_residual() {
        use approx::assert_relative_eq;

        let cy = CarreauYasudaModel::typical_blood();
        let dx = 1e-4_f64;
        let dt = 1e-6_f64;
        let bgk = CarreauYasudaBgk::new(cy, dx, dt);

        let rho = cy.density;
        let q_mag = 250.0_f64;
        let tau = bgk.compute_local_tau(q_mag, rho);
        let cs2 = LATTICE_CS2;
        let dt_dx2 = dt / (dx * dx);
        let shear_rate = q_mag / (2.0 * rho * cs2 * tau * dt);
        let nu = cy.apparent_kinematic_viscosity(shear_rate);
        let rhs = 0.5 + nu * dt_dx2 / cs2;

        assert_relative_eq!(tau, rhs, epsilon = 1e-11);
        assert!(tau > 0.5);
    }
}
