//! Spalart-Allmaras one-equation turbulence model
//!
//! Reference: Spalart, P. R., & Allmaras, S. R. (1994).
//! "A one-equation turbulence model for aerodynamic flows."
//! AIAA Paper 92-0439, 30th Aerospace Sciences Meeting, Reno, NV.
//!
//! Single transport equation for modified turbulent kinematic viscosity ν̃.
//! Popular for aerospace applications due to good boundary layer prediction.

mod helpers;

use self::helpers::{cbrt, wall_distance_field_2d};
use super::constants::{
    EPSILON_MIN, SA_CB1, SA_CB2, SA_CT1, SA_CT2, SA_CT3, SA_CT4, SA_CV1, SA_CV2, SA_CW1, SA_CW2,
    SA_CW3, SA_KAPPA_SQ, SA_SIGMA,
};
use cfd_core::{
    constants::mathematical::numeric::{ONE, ONE_HALF, TWO},
    error::Result,
};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use tracing::instrument;

/// Spalart-Allmaras turbulence model
///
/// Solves transport equation for ν̃ (modified turbulent kinematic viscosity):
/// Dν̃/Dt = Cb1 * S̃ * ν̃ - Cw1 * fw * (ν̃/d)² + (1/σ) * ∇·[(ν + ν̃) ∇ν̃] + Cb2/σ * |∇ν̃|²
///
/// where:
/// - S̃ = modified vorticity magnitude
/// - d = distance to nearest wall
/// - fw = wall destruction function
#[allow(dead_code)]
pub struct SpalartAllmaras<T: RealField + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Model coefficients
    cb1: T,
    cb2: T,
    cw1: T,
    cw2: T,
    cw3: T,
    cv1: T,
    cv2: T,
    ct1: T,
    ct2: T,
    ct3: T,
    ct4: T,
    sigma: T,
    kappa_sq: T,
}

impl<T: RealField + FromPrimitive + Copy> SpalartAllmaras<T> {
    /// Create new Spalart-Allmaras model
    pub fn new(nx: usize, ny: usize) -> Self {
        Self {
            nx,
            ny,
            cb1: T::from_f64(SA_CB1).unwrap_or_else(T::one),
            cb2: T::from_f64(SA_CB2).unwrap_or_else(T::one),
            cw1: T::from_f64(SA_CW1).unwrap_or_else(T::one),
            cw2: T::from_f64(SA_CW2).unwrap_or_else(T::one),
            cw3: T::from_f64(SA_CW3).unwrap_or_else(T::one),
            cv1: T::from_f64(SA_CV1).unwrap_or_else(T::one),
            cv2: T::from_f64(SA_CV2).unwrap_or_else(T::one),
            ct1: T::from_f64(SA_CT1).unwrap_or_else(T::one),
            ct2: T::from_f64(SA_CT2).unwrap_or_else(T::one),
            ct3: T::from_f64(SA_CT3).unwrap_or_else(T::one),
            ct4: T::from_f64(SA_CT4).unwrap_or_else(T::one),
            sigma: T::from_f64(SA_SIGMA).unwrap_or_else(T::one),
            kappa_sq: T::from_f64(SA_KAPPA_SQ).unwrap_or_else(T::one),
        }
    }

    /// Calculate eddy viscosity from modified viscosity
    ///
    /// νt = ν̃ * fv1
    /// where fv1 = χ³/(χ³ + Cv1³) and χ = ν̃/ν
    #[instrument(skip(self))]
    pub fn eddy_viscosity(&self, nu_tilde: T, molecular_viscosity: T) -> T {
        let chi = nu_tilde / molecular_viscosity.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::one));
        let chi_cubed = chi * chi * chi;
        let cv1_cubed = self.cv1 * self.cv1 * self.cv1;
        let fv1 = chi_cubed / (chi_cubed + cv1_cubed);
        nu_tilde * fv1
    }

    /// Calculate vorticity magnitude
    ///
    /// Ω = √(2 * Ωij * Ωij)
    /// where Ωij = 0.5 * (∂ui/∂xj - ∂uj/∂xi) is the rotation tensor
    #[instrument(skip(self))]
    fn vorticity_magnitude(&self, velocity_gradient: &[[T; 2]; 2]) -> T {
        // Rotation tensor Ωij = 0.5 * (dui/dxj - duj/dxi)
        // Ω12 = 0.5 * (du/dy - dv/dx)
        // Ω21 = 0.5 * (dv/dx - du/dy) = -Ω12
        let omega12 = (velocity_gradient[0][1] - velocity_gradient[1][0])
            * T::from_f64(ONE_HALF).unwrap_or_else(T::one);

        // For 2D: Ω = √(2 * (Ω12² + Ω21²)) = √(2 * 2 * Ω12²) = 2|Ω12|
        let two = T::from_f64(TWO).unwrap_or_else(T::one);
        two * omega12.abs()
    }

    /// Calculate modified vorticity S̃
    ///
    /// S̃ = Ω + (ν̃/κ²d²) * fv2
    /// where fv2 = 1 - χ/(1 + χ * fv1)
    #[instrument(skip(self))]
    fn modified_vorticity(
        &self,
        vorticity: T,
        nu_tilde: T,
        molecular_viscosity: T,
        wall_distance: T,
    ) -> T {
        let chi = nu_tilde / molecular_viscosity.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::one));
        let chi_cubed = chi * chi * chi;
        let cv1_cubed = self.cv1 * self.cv1 * self.cv1;
        let fv1 = chi_cubed / (chi_cubed + cv1_cubed);

        // fv2 = 1 - χ/(1 + χ * fv1)
        let fv2 = T::from_f64(ONE).unwrap_or_else(T::one)
            - chi / (T::from_f64(ONE).unwrap_or_else(T::one) + chi * fv1);

        // S̃ = Ω + ν̃/(κ²d²) * fv2
        let d_sq = wall_distance * wall_distance;
        let modification = (nu_tilde / (self.kappa_sq * d_sq.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::one))))
            * fv2;

        vorticity + modification
    }

    /// Calculate production term
    ///
    /// P = Cb1 * S̃ * ν̃
    #[instrument(skip(self))]
    pub fn production(&self, nu_tilde: T, s_tilde: T) -> T {
        self.cb1 * s_tilde * nu_tilde
    }

    /// Calculate wall destruction function fw
    ///
    /// fw = g * [(1 + Cw3⁶)/(g⁶ + Cw3⁶)]^(1/6)
    /// where g = r + Cw2 * (r⁶ - r) and r = ν̃/(S̃ * κ²d²)
    #[instrument(skip(self))]
    fn wall_destruction_function(&self, nu_tilde: T, s_tilde: T, wall_distance: T) -> T {
        let d_sq = wall_distance * wall_distance;
        let denominator = s_tilde * self.kappa_sq * d_sq;

        // Limit r to prevent division by zero
        let r = if denominator.abs() > T::from_f64(EPSILON_MIN).unwrap_or_else(T::one) {
            nu_tilde / denominator
        } else {
            T::from_f64(EPSILON_MIN).unwrap_or_else(T::one)
        };

        // Limit r to reasonable range (prevents numerical issues)
        let r = r.min(T::from_f64(10.0).unwrap_or_else(T::one));

        // g = r + Cw2 * (r⁶ - r)
        let r_sq = r * r;
        let r_6 = r_sq * r_sq * r_sq;
        let g = r + self.cw2 * (r_6 - r);

        // fw = g * [(1 + Cw3⁶)/(g⁶ + Cw3⁶)]^(1/6)
        let cw3_sq = self.cw3 * self.cw3;
        let cw3_6 = cw3_sq * cw3_sq * cw3_sq;
        let g_sq = g * g;
        let g_6 = g_sq * g_sq * g_sq;

        let one = T::from_f64(ONE).unwrap_or_else(T::one);
        let ratio = (one + cw3_6) / (g_6 + cw3_6);

        // Compute sixth root using helper: x^(1/6) = (x^(1/2))^(1/3)
        let sqrt_ratio = ratio.sqrt();
        let cbrt_sqrt = cbrt(sqrt_ratio);

        g * cbrt_sqrt
    }

    /// Calculate destruction term
    ///
    /// D = Cw1 * fw * (ν̃/d)²
    #[instrument(skip(self))]
    pub fn destruction(&self, nu_tilde: T, wall_distance: T, fw: T) -> T {
        let ratio = nu_tilde / wall_distance.max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::one));
        self.cw1 * fw * ratio * ratio
    }

    /// Calculate trip term for transition
    ///
    /// Typically zero for fully turbulent flows, non-zero near laminar-turbulent transition
    /// Ft = Ct1 * gt * exp(-Ct2 * ωt²/Δν² - Ct3 * √(ωt * Δν) * d²)
    #[instrument(skip(self))]
    pub fn trip_term(&self, _nu_tilde: T, _wall_distance: T) -> T {
        // For fully turbulent flows, trip term is zero
        // Can be activated for transition modeling if needed
        T::zero()
    }

    /// Calculate wall distance field
    ///
    /// Simplified 2D wall distance: minimum distance to domain boundaries
    /// In production code, use Eikonal equation solver for complex geometries
    #[instrument(skip(self, dx, dy))]
    pub fn wall_distance_field(&self, dx: T, dy: T) -> Vec<T> {
        wall_distance_field_2d(self.nx, self.ny, dx, dy)
    }

    /// Update modified turbulent viscosity field
    ///
    /// Solves: ∂ν̃/∂t = P - D + (1/σ)∇·[(ν+ν̃)∇ν̃] + (Cb2/σ)|∇ν̃|²
    #[instrument(skip(self, nu_tilde, velocity, molecular_viscosity, dx, dy, dt))]
    pub fn update(
        &self,
        nu_tilde: &mut [T],
        velocity: &[Vector2<T>],
        molecular_viscosity: T,
        dx: T,
        dy: T,
        dt: T,
    ) -> Result<()> {
        // Calculate wall distances
        let wall_distances = self.wall_distance_field(dx, dy);

        // Temporary storage for new values
        let mut nu_new = vec![T::zero(); nu_tilde.len()];

        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                let idx = j * self.nx + i;

                // Compute velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dx);
                let du_dy = (velocity[idx + self.nx].x - velocity[idx - self.nx].x)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dx);
                let dv_dy = (velocity[idx + self.nx].y - velocity[idx - self.nx].y)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dy);

                let velocity_gradient = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate vorticity
                let vorticity = self.vorticity_magnitude(&velocity_gradient);

                // Calculate modified vorticity
                let s_tilde = self.modified_vorticity(
                    vorticity,
                    nu_tilde[idx],
                    molecular_viscosity,
                    wall_distances[idx],
                );

                // Production term
                let production = self.production(nu_tilde[idx], s_tilde);

                // Wall destruction function
                let fw = self.wall_destruction_function(nu_tilde[idx], s_tilde, wall_distances[idx]);

                // Destruction term
                let destruction = self.destruction(nu_tilde[idx], wall_distances[idx], fw);

                // Diffusion term: (1/σ)∇·[(ν+ν̃)∇ν̃]
                let nu_total_e = molecular_viscosity + nu_tilde[idx + 1];
                let nu_total_w = molecular_viscosity + nu_tilde[idx - 1];
                let nu_total_n = molecular_viscosity + nu_tilde[idx + self.nx];
                let nu_total_s = molecular_viscosity + nu_tilde[idx - self.nx];

                let dnu_dx_e = (nu_tilde[idx + 1] - nu_tilde[idx]) / dx;
                let dnu_dx_w = (nu_tilde[idx] - nu_tilde[idx - 1]) / dx;
                let dnu_dy_n = (nu_tilde[idx + self.nx] - nu_tilde[idx]) / dy;
                let dnu_dy_s = (nu_tilde[idx] - nu_tilde[idx - self.nx]) / dy;

                let diffusion_x = (nu_total_e * dnu_dx_e - nu_total_w * dnu_dx_w) / dx;
                let diffusion_y = (nu_total_n * dnu_dy_n - nu_total_s * dnu_dy_s) / dy;
                let diffusion = (diffusion_x + diffusion_y) / self.sigma;

                // Cross-diffusion term: (Cb2/σ)|∇ν̃|²
                let dnu_dx = (nu_tilde[idx + 1] - nu_tilde[idx - 1])
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dx);
                let dnu_dy = (nu_tilde[idx + self.nx] - nu_tilde[idx - self.nx])
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dy);
                let grad_nu_sq = dnu_dx * dnu_dx + dnu_dy * dnu_dy;
                let cross_diffusion = (self.cb2 / self.sigma) * grad_nu_sq;

                // Trip term (zero for fully turbulent)
                let trip = self.trip_term(nu_tilde[idx], wall_distances[idx]);

                // Explicit time integration
                nu_new[idx] = nu_tilde[idx]
                    + dt * (production - destruction + diffusion + cross_diffusion + trip);

                // Ensure positive values
                nu_new[idx] = nu_new[idx].max(T::zero());
            }
        }

        // Update field
        for idx in 0..nu_tilde.len() {
            nu_tilde[idx] = nu_new[idx];
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(nu_tilde);

        Ok(())
    }

    /// Apply boundary conditions
    #[instrument(skip(self, nu_tilde))]
    fn apply_boundary_conditions(&self, nu_tilde: &mut [T]) {
        // Wall boundaries: ν̃ = 0
        for i in 0..self.nx {
            nu_tilde[i] = T::zero(); // Bottom wall
            nu_tilde[(self.ny - 1) * self.nx + i] = T::zero(); // Top wall
        }

        for j in 0..self.ny {
            nu_tilde[j * self.nx] = T::zero(); // Left wall
            nu_tilde[j * self.nx + (self.nx - 1)] = T::zero(); // Right wall
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_spalart_allmaras_creation() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        assert_eq!(model.nx, 10);
        assert_eq!(model.ny, 10);
    }

    #[test]
    fn test_eddy_viscosity_zero_nu_tilde() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_t = model.eddy_viscosity(0.0, 1e-5);
        assert_eq!(nu_t, 0.0);
    }

    #[test]
    fn test_eddy_viscosity_large_nu_tilde() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu = 1e-5;
        let nu_tilde = 1e-3;
        let nu_t = model.eddy_viscosity(nu_tilde, nu);

        // For χ >> 1, fv1 → 1, so νt → ν̃
        assert!(nu_t > 0.9 * nu_tilde);
        assert!(nu_t <= nu_tilde);
    }

    #[test]
    fn test_vorticity_magnitude() {
        let model = SpalartAllmaras::<f64>::new(10, 10);

        // Test case: uniform rotation
        let velocity_gradient = [[0.0, 1.0], [-1.0, 0.0]];
        let vorticity = model.vorticity_magnitude(&velocity_gradient);

        // Ω = |∂v/∂x - ∂u/∂y| = |-1 - 1| = 2
        assert_relative_eq!(vorticity, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_production_term() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_tilde = 1e-4;
        let s_tilde = 100.0;

        let production = model.production(nu_tilde, s_tilde);

        // P = Cb1 * S̃ * ν̃
        let expected = SA_CB1 * s_tilde * nu_tilde;
        assert_relative_eq!(production, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_destruction_term() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let nu_tilde = 1e-4;
        let wall_distance = 0.01;
        let fw = 1.0;

        let destruction = model.destruction(nu_tilde, wall_distance, fw);

        // D = Cw1 * fw * (ν̃/d)²
        let ratio = nu_tilde / wall_distance;
        let expected = SA_CW1 * fw * ratio * ratio;
        assert_relative_eq!(destruction, expected, epsilon = 1e-10);
    }

    #[test]
    fn test_cbrt_function() {
        // Test cube root computation via helper
        use super::helpers::cbrt;
        assert_relative_eq!(cbrt(8.0), 2.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(27.0), 3.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(1.0), 1.0, epsilon = 1e-8);
        assert_eq!(cbrt(0.0), 0.0);
    }

    #[test]
    fn test_wall_distance_field() {
        let model = SpalartAllmaras::<f64>::new(5, 5);
        let distances = model.wall_distance_field(0.1, 0.1);

        assert_eq!(distances.len(), 25);

        // Corner should have minimum distance
        // Center should have maximum distance
        let center_idx = 2 * 5 + 2;
        let corner_idx = 0;

        assert!(distances[center_idx] > distances[corner_idx]);
        assert_relative_eq!(distances[corner_idx], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_boundary_conditions() {
        let model = SpalartAllmaras::<f64>::new(5, 5);
        let mut nu_tilde = vec![1.0; 25];

        model.apply_boundary_conditions(&mut nu_tilde);

        // Check walls are zero
        for i in 0..5 {
            assert_eq!(nu_tilde[i], 0.0); // Bottom
            assert_eq!(nu_tilde[4 * 5 + i], 0.0); // Top
        }
        for j in 0..5 {
            assert_eq!(nu_tilde[j * 5], 0.0); // Left
            assert_eq!(nu_tilde[j * 5 + 4], 0.0); // Right
        }
    }

    #[test]
    fn test_trip_term_zero() {
        let model = SpalartAllmaras::<f64>::new(10, 10);
        let trip = model.trip_term(1e-4, 0.01);

        // For fully turbulent flows, trip term should be zero
        assert_eq!(trip, 0.0);
    }
}
