//! k-ω SST model struct and `TurbulenceModel` trait implementation.
//!
//! # Turbulent Viscosity (Bradshaw Assumption)
//!
//! The SST eddy viscosity is computed with a limiter to ensure the
//! Bradshaw assumption holds in the logarithmic layer:
//!
//! ```text
//! ν_t = a₁ k / max(a₁ ω, S F₂)
//! ```
//!
//! This bounds `τ / ρk ≤ a₁` in the logarithmic and wake regions,
//! preventing overshoot of the Reynolds shear stress.
//!
//! # References
//! - Menter, F.R. (1994). AIAA Journal, 32(8), 1598-1605.
//! - Menter, F.R., Kuntz, M. & Langtry, R. (2003). Turbulence,
//!   Heat and Mass Transfer 4:625-632.

use super::blending::{blend_coefficient, compute_blending_functions, cross_diffusion};
use crate::physics::turbulence::constants::{
    K_MIN, OMEGA_MIN, OMEGA_WALL_COEFFICIENT, SST_ALPHA_1, SST_BETA_1, SST_BETA_2,
    SST_BETA_STAR, SST_GAMMA_1, SST_GAMMA_2, SST_SIGMA_K1, SST_SIGMA_K2, SST_SIGMA_OMEGA1,
    SST_SIGMA_OMEGA2,
};
use crate::physics::turbulence::traits::TurbulenceModel;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// k-ω SST (Shear Stress Transport) turbulence model.
///
/// Combines the robust near-wall k-ω model with the freestream-independent
/// k-ε model through blending functions F1 and F2 (Menter 1994).
pub struct KOmegaSSTModel<T: RealField + Copy> {
    /// Grid x-dimension.
    nx: usize,
    /// Grid y-dimension.
    ny: usize,
    /// Blending function F1 values (flat, `nx * ny`).
    f1: Vec<T>,
    /// Blending function F2 values (flat, `nx * ny`).
    f2: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> KOmegaSSTModel<T> {
    /// Create a new k-ω SST model for a grid of dimension `nx × ny`.
    #[must_use]
    pub fn new(nx: usize, ny: usize) -> Self {
        let size = nx * ny;
        Self {
            nx,
            ny,
            f1: vec![T::zero(); size],
            f2: vec![T::zero(); size],
        }
    }

    /// Apply wall and positivity boundary conditions.
    fn apply_boundary_conditions(&self, k: &mut [T], omega: &mut [T], wall_distance: &[T]) {
        let omega_min = T::from_f64(OMEGA_MIN).expect("analytical constant conversion");

        for i in 0..k.len() {
            k[i] = k[i].max(omega_min);
            omega[i] = omega[i].max(omega_min);
        }

        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;
                const WALL_PROXIMITY_THRESHOLD: f64 = 1e-6;
                if wall_distance[idx]
                    < T::from_f64(WALL_PROXIMITY_THRESHOLD).expect("analytical constant conversion")
                {
                    k[idx] = T::zero();
                    let y = wall_distance[idx].max(T::from_f64(1e-10).expect("analytical constant conversion"));
                    omega[idx] =
                        T::from_f64(OMEGA_WALL_COEFFICIENT).expect("analytical constant conversion") / (y * y);
                }
            }
        }
    }
}

impl<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive> TurbulenceModel<T>
    for KOmegaSSTModel<T>
{
    fn turbulent_viscosity(&self, k: T, omega: T, density: T) -> T {
        let omega_min = T::from_f64(OMEGA_MIN).expect("analytical constant conversion");
        let a1 = T::from_f64(SST_ALPHA_1).expect("analytical constant conversion");
        let nu_t_unlimited = k / omega.max(omega_min);
        density * nu_t_unlimited.min(a1 * k / omega.max(omega_min))
    }

    fn turbulent_viscosity_with_limiter(
        &self,
        k: T,
        omega: T,
        density: T,
        strain_rate_magnitude: T,
        f2: T,
    ) -> T {
        let omega_min = T::from_f64(OMEGA_MIN).expect("analytical constant conversion");
        let a1 = T::from_f64(SST_ALPHA_1).expect("analytical constant conversion");
        let denominator = (a1 * omega).max(strain_rate_magnitude * f2);
        let nu_t = a1 * k / denominator.max(omega_min);
        density * nu_t
    }

    fn production_term(
        &self,
        velocity_gradient: &[[T; 2]; 2],
        turbulent_viscosity: T,
        _turbulence_variable: T,
        _wall_distance: T,
        _molecular_viscosity: T,
    ) -> T {
        let mut s_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                let s_ij = (velocity_gradient[i][j] + velocity_gradient[j][i])
                    * T::from_f64(0.5).expect("analytical constant conversion");
                s_squared += s_ij * s_ij;
            }
        }
        let two = T::from_f64(2.0).expect("analytical constant conversion");
        turbulent_viscosity * two * s_squared
    }

    fn dissipation_term(&self, k: T, omega: T) -> T {
        let beta_star = T::from_f64(SST_BETA_STAR).expect("analytical constant conversion");
        beta_star * k * omega
    }

    fn update(
        &mut self,
        k: &mut [T],
        omega: &mut [T],
        velocity: &[Vector2<T>],
        density: T,
        molecular_viscosity: T,
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;

        // Wall distance (channel geometry: walls at y=0 and y=H)
        let mut wall_distance = vec![T::one(); nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                let y_bottom = T::from_usize(j).expect("analytical constant conversion") * dy;
                let y_top = T::from_usize(ny - 1 - j).expect("analytical constant conversion") * dy;
                wall_distance[idx] = y_bottom.min(y_top);
            }
        }

        // Compute blending functions
        compute_blending_functions(
            &mut self.f1,
            &mut self.f2,
            k,
            omega,
            &wall_distance,
            molecular_viscosity,
            density,
            nx,
            ny,
            dx,
            dy,
        );

        let k_previous = k.to_vec();
        let omega_previous = omega.to_vec();

        // Update interior points
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                let beta = blend_coefficient(SST_BETA_1, SST_BETA_2, self.f1[idx]);
                let gamma = blend_coefficient(SST_GAMMA_1, SST_GAMMA_2, self.f1[idx]);
                let sigma_k = blend_coefficient(SST_SIGMA_K1, SST_SIGMA_K2, self.f1[idx]);
                let sigma_omega =
                    blend_coefficient(SST_SIGMA_OMEGA1, SST_SIGMA_OMEGA2, self.f1[idx]);

                // Velocity gradients (central differences)
                let two = T::from_f64(2.0).expect("analytical constant conversion");
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x) / (two * dx);
                let du_dy = (velocity[idx + nx].x - velocity[idx - nx].x) / (two * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y) / (two * dx);
                let dv_dy = (velocity[idx + nx].y - velocity[idx - nx].y) / (two * dy);
                let grad = [[du_dx, du_dy], [dv_dx, dv_dy]];

                let nu_t =
                    self.turbulent_viscosity(k_previous[idx], omega_previous[idx], density);

                // Production with Menter (2003) limiter
                let p_k_unlimited = self.production_term(
                    &grad,
                    nu_t,
                    k_previous[idx],
                    wall_distance[idx],
                    molecular_viscosity,
                );
                let c_lim = T::from_f64(10.0).expect("analytical constant conversion");
                let beta_star = T::from_f64(SST_BETA_STAR).expect("analytical constant conversion");
                let p_k = p_k_unlimited
                    .min(c_lim * beta_star * k_previous[idx] * omega_previous[idx]);

                // Effective viscosities
                let nu_eff_k = molecular_viscosity + nu_t * sigma_k;
                let nu_eff_omega = molecular_viscosity + nu_t * sigma_omega;

                // k diffusion (2nd-order central)
                let diff_k_x = (k_previous[idx + 1] - two * k_previous[idx]
                    + k_previous[idx - 1])
                    / (dx * dx);
                let diff_k_y = (k_previous[idx + nx] - two * k_previous[idx]
                    + k_previous[idx - nx])
                    / (dy * dy);
                let diff_k = nu_eff_k * (diff_k_x + diff_k_y);

                // ω diffusion
                let diff_omega_x = (omega_previous[idx + 1] - two * omega_previous[idx]
                    + omega_previous[idx - 1])
                    / (dx * dx);
                let diff_omega_y = (omega_previous[idx + nx] - two * omega_previous[idx]
                    + omega_previous[idx - nx])
                    / (dy * dy);
                let diff_omega = nu_eff_omega * (diff_omega_x + diff_omega_y);

                // Cross-diffusion
                let cd_kw =
                    cross_diffusion(&k_previous, &omega_previous, idx, nx, ny, dx, dy);
                let cd_term = (T::one() - self.f1[idx]) * cd_kw;

                // Advance k
                let k_new = k_previous[idx]
                    + dt * (p_k - beta_star * k_previous[idx] * omega_previous[idx] + diff_k);
                let k_min = T::from_f64(K_MIN).expect("analytical constant conversion");
                k[idx] = k_new.max(k_min);

                // Advance ω
                let omega_source = gamma * density * p_k
                    / nu_t.max(T::from_f64(OMEGA_MIN).expect("analytical constant conversion"));
                let omega_sink = beta * density * omega_previous[idx] * omega_previous[idx];
                let omega_new =
                    omega_previous[idx] + dt * (omega_source - omega_sink + diff_omega + cd_term);
                let omega_min = T::from_f64(OMEGA_MIN).expect("analytical constant conversion");
                omega[idx] = omega_new.max(omega_min);
            }
        }

        self.apply_boundary_conditions(k, omega, &wall_distance);
        Ok(())
    }

    fn name(&self) -> &'static str {
        "k-omega-SST"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        reynolds > T::zero()
    }
}
