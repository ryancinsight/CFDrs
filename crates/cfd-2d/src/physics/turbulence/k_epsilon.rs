//! k-ε turbulence model implementation
//!
//! Based on Launder & Spalding (1974) "The numerical computation of turbulent flows"
//! Computer Methods in Applied Mechanics and Engineering, 3(2), 269-289
//!
//! Standard k-ε model with wall functions for industrial CFD applications

use super::constants::{
    C1_EPSILON, C2_EPSILON, C_MU, EPSILON_MIN, SIGMA_EPSILON, SIGMA_K,
};
use super::traits::TurbulenceModel;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;

/// k-ε turbulence model
pub struct KEpsilonModel<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Model coefficients (can be modified for realizability)
    c_mu: T,
    c1_epsilon: T,
    c2_epsilon: T,
    sigma_k: T,
    sigma_epsilon: T,
}

impl<T: RealField + FromPrimitive + Copy> KEpsilonModel<T> {
    /// Create a new k-ε model
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            c_mu: T::from_f64(C_MU).unwrap_or_else(T::one),
            c1_epsilon: T::from_f64(C1_EPSILON).unwrap_or_else(T::one),
            c2_epsilon: T::from_f64(C2_EPSILON).unwrap_or_else(T::one),
            sigma_k: T::from_f64(SIGMA_K).unwrap_or_else(T::one),
            sigma_epsilon: T::from_f64(SIGMA_EPSILON).unwrap_or_else(T::one),
        }
    }

    /// Calculate strain rate tensor
    fn strain_rate(&self, velocity_gradient: &[[T; 2]; 2]) -> T {
        let mut s_ij = [[T::zero(); 2]; 2];

        // S_ij = 0.5 * (du_i/dx_j + du_j/dx_i)
        for i in 0..2 {
            for j in 0..2 {
                s_ij[i][j] = (velocity_gradient[i][j] + velocity_gradient[j][i])
                    * T::from_f64(super::constants::HALF).unwrap_or_else(T::one);
            }
        }

        // Calculate magnitude: sqrt(2 * S_ij * S_ij)
        let mut s_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                s_squared += s_ij[i][j] * s_ij[i][j];
            }
        }

        (T::from_f64(super::constants::TWO).unwrap_or_else(T::one) * s_squared).sqrt()
    }

    /// Apply boundary conditions
    fn apply_boundary_conditions(&self, k: &mut [T], epsilon: &mut [T]) {
        let eps_min = T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero);

        // Ensure positive values
        for i in 0..k.len() {
            k[i] = k[i].max(eps_min);
            epsilon[i] = epsilon[i].max(eps_min);
        }

        // Wall boundaries (assuming walls at j=0 and j=ny-1)
        for i in 0..self.nx {
            // Bottom wall
            k[i] = T::zero();
            epsilon[i] = eps_min;

            // Top wall
            let idx_top = i + (self.ny - 1) * self.nx;
            k[idx_top] = T::zero();
            epsilon[idx_top] = eps_min;
        }
    }
}

impl<T: RealField + FromPrimitive + Copy> TurbulenceModel<T> for KEpsilonModel<T> {
    fn turbulent_viscosity(&self, k: T, epsilon: T, density: T) -> T {
        let eps_min = T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero);
        density * self.c_mu * k * k / epsilon.max(eps_min)
    }

    fn production_term(&self, velocity_gradient: &[[T; 2]; 2], turbulent_viscosity: T) -> T {
        let strain = self.strain_rate(velocity_gradient);
        let two = T::from_f64(super::constants::TWO).unwrap_or_else(T::one);
        turbulent_viscosity * strain * strain * two
    }

    fn dissipation_term(&self, k: T, epsilon: T) -> T {
        epsilon
    }

    fn update(
        &mut self,
        k: &mut [T],
        epsilon: &mut [T],
        velocity: &[Vector2<T>],
        density: T,
        molecular_viscosity: T,
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;

        // Store old values for explicit time stepping
        let k_old = k.to_vec();
        let epsilon_old = epsilon.to_vec();

        // Update interior points
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                // Calculate velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x)
                    / (T::from_f64(super::constants::TWO).unwrap_or_else(T::one) * dx);
                let du_dy = (velocity[idx + nx].x - velocity[idx - nx].x)
                    / (T::from_f64(super::constants::TWO).unwrap_or_else(T::one) * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y)
                    / (T::from_f64(super::constants::TWO).unwrap_or_else(T::one) * dx);
                let dv_dy = (velocity[idx + nx].y - velocity[idx - nx].y)
                    / (T::from_f64(super::constants::TWO).unwrap_or_else(T::one) * dy);

                let grad = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate turbulent viscosity
                let nu_t = self.turbulent_viscosity(k_old[idx], epsilon_old[idx], density);

                // Production term
                let p_k = self.production_term(&grad, nu_t);

                // Diffusion terms
                let nu_eff_k = molecular_viscosity + nu_t / self.sigma_k;
                let nu_eff_eps = molecular_viscosity + nu_t / self.sigma_epsilon;

                // k equation diffusion
                let diff_k_x = (k_old[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_old[idx]
                    + k_old[idx - 1])
                    / (dx * dx);
                let diff_k_y = (k_old[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_old[idx]
                    + k_old[idx - nx])
                    / (dy * dy);
                let diff_k = nu_eff_k * (diff_k_x + diff_k_y);

                // epsilon equation diffusion
                let diff_eps_x = (epsilon_old[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * epsilon_old[idx]
                    + epsilon_old[idx - 1])
                    / (dx * dx);
                let diff_eps_y = (epsilon_old[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * epsilon_old[idx]
                    + epsilon_old[idx - nx])
                    / (dy * dy);
                let diff_eps = nu_eff_eps * (diff_eps_x + diff_eps_y);

                // Update k
                k[idx] = k_old[idx] + dt * (p_k - epsilon_old[idx] + diff_k);

                // Update epsilon
                let eps_source = self.c1_epsilon * epsilon_old[idx]
                    / k_old[idx].max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero))
                    * p_k;
                let eps_sink = self.c2_epsilon * epsilon_old[idx] * epsilon_old[idx]
                    / k_old[idx].max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero));
                epsilon[idx] = epsilon_old[idx] + dt * (eps_source - eps_sink + diff_eps);
            }
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(k, epsilon);

        Ok(())
    }

    fn name(&self) -> &str {
        "k-epsilon"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        // k-ε is valid for high Reynolds numbers
        reynolds > T::from_f64(1e4).unwrap_or_else(T::one)
    }
}
