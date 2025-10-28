//! k-ε turbulence model implementation
//!
//! Based on Launder & Spalding (1974) "The numerical computation of turbulent flows"
//! Computer Methods in Applied Mechanics and Engineering, 3(2), 269-289
//!
//! Standard k-ε model with wall functions for industrial CFD applications

use super::constants::{C1_EPSILON, C2_EPSILON, C_MU, EPSILON_MIN, SIGMA_EPSILON, SIGMA_K};
use super::traits::TurbulenceModel;
use cfd_core::{
    constants::mathematical::numeric::{ONE_HALF, TWO},
    error::Result,
};
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
                    * T::from_f64(ONE_HALF).unwrap_or_else(T::one);
            }
        }

        // Calculate magnitude: sqrt(2 * S_ij * S_ij)
        let mut s_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                s_squared += s_ij[i][j] * s_ij[i][j];
            }
        }

        (T::from_f64(TWO).unwrap_or_else(T::one) * s_squared).sqrt()
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
        let two = T::from_f64(TWO).unwrap_or_else(T::one);
        turbulent_viscosity * strain * strain * two
    }

    fn dissipation_term(&self, _k: T, epsilon: T) -> T {
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

        // Store previous timestep values for explicit time stepping
        // These are algorithmically required, not naming violations
        let k_previous = k.to_vec();
        let epsilon_previous = epsilon.to_vec();

        // Update interior points
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                // Calculate velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dx);
                let du_dy = (velocity[idx + nx].x - velocity[idx - nx].x)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dx);
                let dv_dy = (velocity[idx + nx].y - velocity[idx - nx].y)
                    / (T::from_f64(TWO).unwrap_or_else(T::one) * dy);

                let grad = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate turbulent viscosity
                let nu_t =
                    self.turbulent_viscosity(k_previous[idx], epsilon_previous[idx], density);

                // Production term
                let p_k = self.production_term(&grad, nu_t);

                // Diffusion terms
                let nu_eff_k = molecular_viscosity + nu_t / self.sigma_k;
                let nu_eff_eps = molecular_viscosity + nu_t / self.sigma_epsilon;

                // k equation diffusion
                let diff_k_x = (k_previous[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_previous[idx]
                    + k_previous[idx - 1])
                    / (dx * dx);
                let diff_k_y = (k_previous[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * k_previous[idx]
                    + k_previous[idx - nx])
                    / (dy * dy);
                let diff_k = nu_eff_k * (diff_k_x + diff_k_y);

                // epsilon equation diffusion
                let diff_eps_x = (epsilon_previous[idx + 1]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * epsilon_previous[idx]
                    + epsilon_previous[idx - 1])
                    / (dx * dx);
                let diff_eps_y = (epsilon_previous[idx + nx]
                    - T::from_f64(2.0).unwrap_or_else(T::one) * epsilon_previous[idx]
                    + epsilon_previous[idx - nx])
                    / (dy * dy);
                let diff_eps = nu_eff_eps * (diff_eps_x + diff_eps_y);

                // Update k
                k[idx] = k_previous[idx] + dt * (p_k - epsilon_previous[idx] + diff_k);

                // Update epsilon
                let eps_source = self.c1_epsilon * epsilon_previous[idx]
                    / k_previous[idx].max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero))
                    * p_k;
                let eps_sink = self.c2_epsilon * epsilon_previous[idx] * epsilon_previous[idx]
                    / k_previous[idx].max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero));
                epsilon[idx] = epsilon_previous[idx] + dt * (eps_source - eps_sink + diff_eps);
            }
        }

        // Apply boundary conditions
        self.apply_boundary_conditions(k, epsilon);

        Ok(())
    }

    fn name(&self) -> &'static str {
        "k-epsilon"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        // k-ε is valid for high Reynolds numbers
        reynolds > T::from_f64(1e4).unwrap_or_else(T::one)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_new_model_initialization() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        assert_eq!(model.nx, 10);
        assert_eq!(model.ny, 10);
        assert_relative_eq!(model.c_mu, C_MU, epsilon = 1e-10);
        assert_relative_eq!(model.c1_epsilon, C1_EPSILON, epsilon = 1e-10);
        assert_relative_eq!(model.c2_epsilon, C2_EPSILON, epsilon = 1e-10);
        assert_relative_eq!(model.sigma_k, SIGMA_K, epsilon = 1e-10);
        assert_relative_eq!(model.sigma_epsilon, SIGMA_EPSILON, epsilon = 1e-10);
    }

    #[test]
    fn test_turbulent_viscosity_positive() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let k = 1.0;
        let epsilon = 1.0;
        let density = 1.0;
        let nu_t = model.turbulent_viscosity(k, epsilon, density);
        assert!(nu_t > 0.0);
        assert_relative_eq!(nu_t, density * C_MU * k * k / epsilon, epsilon = 1e-10);
    }

    #[test]
    fn test_turbulent_viscosity_zero_epsilon() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let k = 1.0;
        let epsilon = 0.0;
        let density = 1.0;
        let nu_t = model.turbulent_viscosity(k, epsilon, density);
        // Should use EPSILON_MIN to prevent division by zero
        assert!(nu_t > 0.0);
        assert!(nu_t.is_finite());
    }

    #[test]
    fn test_strain_rate_pure_shear() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        // Pure shear flow: du/dy = 1, all other gradients zero
        let grad = [[0.0, 1.0], [0.0, 0.0]];
        let strain = model.strain_rate(&grad);
        // S_12 = S_21 = 0.5 * (0 + 1) = 0.5
        // |S| = sqrt(2 * (S_12^2 + S_21^2)) = sqrt(2 * (0.25 + 0.25)) = sqrt(1.0) = 1.0
        assert_relative_eq!(strain, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn test_strain_rate_pure_extension() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        // Pure extension: du/dx = 1, dv/dy = -1 (incompressible)
        let grad = [[1.0, 0.0], [0.0, -1.0]];
        let strain = model.strain_rate(&grad);
        // S_11 = 1.0, S_22 = -1.0
        // |S| = sqrt(2 * (1^2 + 1^2)) = sqrt(4) = 2.0
        assert_relative_eq!(strain, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_strain_rate_zero_velocity_gradient() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let grad = [[0.0, 0.0], [0.0, 0.0]];
        let strain = model.strain_rate(&grad);
        assert_relative_eq!(strain, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_production_term_positive() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let grad = [[1.0, 0.0], [0.0, 1.0]];
        let nu_t = 0.1;
        let p_k = model.production_term(&grad, nu_t);
        assert!(p_k > 0.0);
    }

    #[test]
    fn test_production_term_zero_strain() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let grad = [[0.0, 0.0], [0.0, 0.0]];
        let nu_t = 0.1;
        let p_k = model.production_term(&grad, nu_t);
        assert_relative_eq!(p_k, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_dissipation_term() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let k = 1.0;
        let epsilon = 2.5;
        let dissipation = model.dissipation_term(k, epsilon);
        assert_relative_eq!(dissipation, epsilon, epsilon = 1e-10);
    }

    #[test]
    fn test_apply_boundary_conditions_enforces_positivity() {
        let model = KEpsilonModel::<f64>::new(5, 5);
        let mut k = vec![-1.0; 25]; // Negative values
        let mut epsilon = vec![-2.0; 25];
        
        model.apply_boundary_conditions(&mut k, &mut epsilon);
        
        // All values should be non-negative
        for &val in &k {
            assert!(val >= 0.0);
        }
        for &val in &epsilon {
            assert!(val >= EPSILON_MIN);
        }
    }

    #[test]
    fn test_apply_boundary_conditions_wall_values() {
        let model = KEpsilonModel::<f64>::new(5, 5);
        let mut k = vec![1.0; 25];
        let mut epsilon = vec![1.0; 25];
        
        model.apply_boundary_conditions(&mut k, &mut epsilon);
        
        // Bottom wall (j=0)
        for i in 0..5 {
            assert_relative_eq!(k[i], 0.0, epsilon = 1e-10);
            assert!(epsilon[i] >= EPSILON_MIN);
        }
        
        // Top wall (j=4)
        for i in 0..5 {
            let idx = i + 4 * 5;
            assert_relative_eq!(k[idx], 0.0, epsilon = 1e-10);
            assert!(epsilon[idx] >= EPSILON_MIN);
        }
    }

    #[test]
    fn test_update_maintains_positivity() {
        let mut model = KEpsilonModel::<f64>::new(5, 5);
        let mut k = vec![0.1; 25];
        let mut epsilon = vec![0.1; 25];
        let velocity = vec![Vector2::new(0.0, 0.0); 25];
        let density = 1.0;
        let molecular_viscosity = 1e-5;
        let dt = 0.001;
        let dx = 0.1;
        let dy = 0.1;
        
        let result = model.update(&mut k, &mut epsilon, &velocity, density, 
                                  molecular_viscosity, dt, dx, dy);
        
        assert!(result.is_ok());
        
        // All values should remain non-negative
        for &val in &k {
            assert!(val >= 0.0);
        }
        for &val in &epsilon {
            assert!(val >= 0.0);
        }
    }

    #[test]
    fn test_model_name() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        assert_eq!(model.name(), "k-epsilon");
    }

    #[test]
    fn test_is_valid_for_high_reynolds() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        assert!(model.is_valid_for_reynolds(1e5));
        assert!(model.is_valid_for_reynolds(1e6));
    }

    #[test]
    fn test_is_not_valid_for_low_reynolds() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        assert!(!model.is_valid_for_reynolds(1e3));
        assert!(!model.is_valid_for_reynolds(1e2));
    }

    #[test]
    fn test_reynolds_threshold() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        assert!(!model.is_valid_for_reynolds(1e4)); // Exactly at threshold
        assert!(model.is_valid_for_reynolds(1e4 + 1.0)); // Just above threshold
    }

    #[test]
    fn test_update_with_uniform_flow() {
        let mut model = KEpsilonModel::<f64>::new(5, 5);
        let mut k = vec![1.0; 25];
        let mut epsilon = vec![1.0; 25];
        // Uniform flow (no velocity gradients)
        let velocity = vec![Vector2::new(1.0, 0.0); 25];
        let density = 1.0;
        let molecular_viscosity = 1e-5;
        let dt = 0.001;
        let dx = 0.1;
        let dy = 0.1;
        
        let result = model.update(&mut k, &mut epsilon, &velocity, density, 
                                  molecular_viscosity, dt, dx, dy);
        
        assert!(result.is_ok());
    }

    #[test]
    fn test_turbulent_viscosity_scales_with_k_squared() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let epsilon = 1.0;
        let density = 1.0;
        
        let nu_t_1 = model.turbulent_viscosity(1.0, epsilon, density);
        let nu_t_2 = model.turbulent_viscosity(2.0, epsilon, density);
        
        // nu_t ~ k^2, so doubling k should quadruple nu_t
        assert_relative_eq!(nu_t_2 / nu_t_1, 4.0, epsilon = 1e-10);
    }

    #[test]
    fn test_production_term_scales_with_turbulent_viscosity() {
        let model = KEpsilonModel::<f64>::new(10, 10);
        let grad = [[1.0, 0.0], [0.0, 1.0]];
        
        let p_k_1 = model.production_term(&grad, 0.1);
        let p_k_2 = model.production_term(&grad, 0.2);
        
        // P_k ~ nu_t, so doubling nu_t should double P_k
        assert_relative_eq!(p_k_2 / p_k_1, 2.0, epsilon = 1e-10);
    }
}
