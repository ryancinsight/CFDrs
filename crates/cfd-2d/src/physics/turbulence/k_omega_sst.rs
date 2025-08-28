//! k-ω SST turbulence model implementation

use super::constants::*;
use super::traits::TurbulenceModel;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// k-ω SST (Shear Stress Transport) turbulence model
pub struct KOmegaSSTModel<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Blending function values
    f1: Vec<T>,
    f2: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> KOmegaSSTModel<T> {
    /// Create a new k-ω SST model
    pub fn new(nx: usize, ny: usize) -> Self {
        let size = nx * ny;
        Self {
            nx,
            ny,
            f1: vec![T::zero(); size],
            f2: vec![T::zero(); size],
        }
    }

    /// Calculate cross-diffusion term CDkω
    fn calculate_cross_diffusion(&self, k: &[T], omega: &[T], idx: usize, dx: T, dy: T) -> T {
        let nx = self.nx;
        let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);

        // Calculate gradients using central differences
        let i = idx % nx;
        let j = idx / nx;

        // Boundary check
        if i == 0 || i == nx - 1 || j == 0 || j == self.ny - 1 {
            return omega_min;
        }

        // ∇k
        let dk_dx = (k[idx + 1] - k[idx - 1]) / (T::from_f64(2.0).unwrap_or_else(T::one) * dx);
        let dk_dy = (k[idx + nx] - k[idx - nx]) / (T::from_f64(2.0).unwrap_or_else(T::one) * dy);

        // ∇ω
        let domega_dx =
            (omega[idx + 1] - omega[idx - 1]) / (T::from_f64(2.0).unwrap_or_else(T::one) * dx);
        let domega_dy =
            (omega[idx + nx] - omega[idx - nx]) / (T::from_f64(2.0).unwrap_or_else(T::one) * dy);

        // Dot product ∇k · ∇ω
        let grad_dot = dk_dx * domega_dx + dk_dy * domega_dy;

        // CDkω = 2ρσω2 / ω · ∇k · ∇ω
        let sigma_omega2 = T::from_f64(SST_SIGMA_OMEGA2).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);

        (two * sigma_omega2 * grad_dot / omega[idx].max(omega_min)).max(omega_min)
    }

    /// Calculate blending functions F1 and F2
    fn calculate_blending_functions(
        &mut self,
        k: &[T],
        omega: &[T],
        wall_distance: &[T],
        molecular_viscosity: T,
        density: T,
        dx: T,
        dy: T,
    ) {
        let beta_star = T::from_f64(SST_BETA_STAR).unwrap_or_else(T::one);
        let sigma_omega2 = T::from_f64(SST_SIGMA_OMEGA2).unwrap_or_else(T::one);
        let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);

        for idx in 0..k.len() {
            let y = wall_distance[idx];
            let k_val = k[idx].max(omega_min);
            let omega_val = omega[idx].max(omega_min);

            // Calculate vorticity magnitude for 2D
            let nu = molecular_viscosity / density;

            // CDkω calculation (cross-diffusion term)
            let cd_kw = self.calculate_cross_diffusion(k, omega, idx, dx, dy);

            // Argument for F1
            let sqrt_k = k_val.sqrt();
            let arg1_1 = sqrt_k / (beta_star * omega_val * y);
            let arg1_2 = T::from_f64(500.0).unwrap_or_else(T::one) * nu / (y * y * omega_val);
            let arg1_3 = T::from_f64(4.0).unwrap_or_else(T::one) * density * sigma_omega2 * k_val
                / (cd_kw * y * y);

            let arg1 = arg1_1.min(arg1_2).max(arg1_3);
            self.f1[idx] = (T::from_f64(4.0).unwrap_or_else(T::one) * arg1).tanh();

            // Argument for F2
            let arg2_1 =
                T::from_f64(2.0).unwrap_or_else(T::one) * sqrt_k / (beta_star * omega_val * y);
            let arg2_2 = T::from_f64(500.0).unwrap_or_else(T::one) * nu / (y * y * omega_val);

            let arg2 = arg2_1.max(arg2_2);
            self.f2[idx] = (arg2 * arg2).tanh();
        }
    }

    /// Blend coefficients based on F1
    fn blend_coefficient(&self, coef1: f64, coef2: f64, f1: T) -> T {
        let c1 = T::from_f64(coef1).unwrap_or_else(T::one);
        let c2 = T::from_f64(coef2).unwrap_or_else(T::one);
        f1 * c1 + (T::one() - f1) * c2
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&self, k: &mut [T], omega: &mut [T], wall_distance: &[T]) {
        let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);

        // Ensure positive values
        for i in 0..k.len() {
            k[i] = k[i].max(omega_min);
            omega[i] = omega[i].max(omega_min);
        }

        // Wall boundaries with omega wall treatment
        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;

                // Near-wall treatment based on Menter (1994) SST model
                // Wall proximity threshold from DNS studies
                const WALL_PROXIMITY_THRESHOLD: f64 = 1e-6;
                if wall_distance[idx]
                    < T::from_f64(WALL_PROXIMITY_THRESHOLD).unwrap_or_else(T::zero)
                {
                    k[idx] = T::zero();
                    // Omega wall value: 60*nu/(beta*y^2)
                    let y = wall_distance[idx].max(T::from_f64(1e-10).unwrap_or_else(T::zero));
                    omega[idx] =
                        T::from_f64(OMEGA_WALL_COEFFICIENT).unwrap_or_else(T::one) / (y * y);
                }
            }
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> TurbulenceModel<T> for KOmegaSSTModel<T> {
    fn turbulent_viscosity(&self, k: T, omega: T, density: T) -> T {
        let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);
        let a1 = T::from_f64(SST_ALPHA_1).unwrap_or_else(T::one);

        // SST limiter for eddy viscosity
        let nu_t_unlimited = k / omega.max(omega_min);

        // Apply SST limiter - Bradshaw assumption
        // νt = a1*k / max(a1*ω, S*F2)
        // For now using standard limiter without strain rate
        // Full implementation requires strain rate magnitude S
        density * nu_t_unlimited.min(a1 * k / omega.max(omega_min))
    }

    fn production_term(&self, velocity_gradient: &[[T; 2]; 2], turbulent_viscosity: T) -> T {
        // Calculate strain rate tensor magnitude
        let mut s_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                let s_ij = (velocity_gradient[i][j] + velocity_gradient[j][i])
                    * T::from_f64(0.5).unwrap_or_else(T::one);
                s_squared += s_ij * s_ij;
            }
        }

        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        turbulent_viscosity * two * s_squared
    }

    fn dissipation_term(&self, k: T, omega: T) -> T {
        let beta_star = T::from_f64(SST_BETA_STAR).unwrap_or_else(T::one);
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

        // Calculate wall distance
        // Using geometric distance for channel flow (walls at y=0 and y=H)
        // For complex geometries, use Poisson equation or fast marching method
        let mut wall_distance = vec![T::one(); nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;
                // Distance from bottom wall
                let y_bottom = T::from_usize(j).unwrap_or_else(T::zero) * dy;
                // Distance from top wall
                let y_top = T::from_usize(ny - 1 - j).unwrap_or_else(T::zero) * dy;
                wall_distance[idx] = y_bottom.min(y_top);
            }
        }

        // Calculate blending functions
        self.calculate_blending_functions(
            k,
            omega,
            &wall_distance,
            molecular_viscosity,
            density,
            dx,
            dy,
        );

        // Store old values
        let k_old = k.to_vec();
        let omega_old = omega.to_vec();

        // Update interior points
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                // Get blended coefficients
                let beta = self.blend_coefficient(SST_BETA_1, SST_BETA_2, self.f1[idx]);
                let gamma = self.blend_coefficient(SST_GAMMA_1, SST_GAMMA_2, self.f1[idx]);
                let sigma_k = self.blend_coefficient(SST_SIGMA_K1, SST_SIGMA_K2, self.f1[idx]);
                let sigma_omega =
                    self.blend_coefficient(SST_SIGMA_OMEGA1, SST_SIGMA_OMEGA2, self.f1[idx]);

                // Calculate velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x)
                    / (T::from_f64(2.0).unwrap_or_else(T::one) * dx);
                let du_dy = (velocity[idx + nx].x - velocity[idx - nx].x)
                    / (T::from_f64(2.0).unwrap_or_else(T::one) * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y)
                    / (T::from_f64(2.0).unwrap_or_else(T::one) * dx);
                let dv_dy = (velocity[idx + nx].y - velocity[idx - nx].y)
                    / (T::from_f64(2.0).unwrap_or_else(T::one) * dy);

                let grad = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate turbulent viscosity
                let nu_t = self.turbulent_viscosity(k_old[idx], omega_old[idx], density);

                // Production term
                let p_k = self.production_term(&grad, nu_t);

                // Diffusion terms
                let nu_eff_k = molecular_viscosity + nu_t * sigma_k;
                let nu_eff_omega = molecular_viscosity + nu_t * sigma_omega;

                // k equation
                let diff_k_x = (k_old[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_old[idx]
                    + k_old[idx - 1])
                    / (dx * dx);
                let diff_k_y = (k_old[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_old[idx]
                    + k_old[idx - nx])
                    / (dy * dy);
                let diff_k = nu_eff_k * (diff_k_x + diff_k_y);

                // omega equation
                let diff_omega_x = (omega_old[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * omega_old[idx]
                    + omega_old[idx - 1])
                    / (dx * dx);
                let diff_omega_y = (omega_old[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * omega_old[idx]
                    + omega_old[idx - nx])
                    / (dy * dy);
                let diff_omega = nu_eff_omega * (diff_omega_x + diff_omega_y);

                // Cross-diffusion term
                let cd_kw = self.calculate_cross_diffusion(&k_old, &omega_old, idx, dx, dy);
                let cd_term = (T::one() - self.f1[idx]) * cd_kw;

                // Update k
                let beta_star = T::from_f64(SST_BETA_STAR).unwrap_or_else(T::one);
                k[idx] = k_old[idx] + dt * (p_k - beta_star * k_old[idx] * omega_old[idx] + diff_k);

                // Update omega
                let omega_source = gamma * density * p_k
                    / nu_t.max(T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero));
                let omega_sink = beta * density * omega_old[idx] * omega_old[idx];
                omega[idx] =
                    omega_old[idx] + dt * (omega_source - omega_sink + diff_omega + cd_term);
            }
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(k, omega, &wall_distance);

        Ok(())
    }

    fn name(&self) -> &str {
        "k-omega-SST"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        // SST is valid for all Reynolds numbers
        reynolds > T::zero()
    }
}
