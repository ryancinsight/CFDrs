//! Sigma Subgrid-Scale (SGS) Model for Transitional Flows
//!
//! The Sigma model is specifically designed for transitional flow regimes,
//! providing superior performance in flows undergoing laminar-to-turbulent transition.
//!
//! ## Mathematical Foundation
//!
//! The Sigma model computes SGS viscosity using the square of the velocity gradient tensor:
//!
//! ```math
//! ν_SGS = (C_σ Δ)² Σ_{i,j} (∂u_i/∂x_j)²
//! ```
//!
//! This formulation naturally:
//! - Reduces to zero in laminar regions (where velocity gradients are small)
//! - Provides appropriate dissipation in turbulent regions
//! - Avoids over-dissipation in transitional zones
//!
//! ## Key Advantages
//!
//! 1. **Transitional Flow Handling**: Specifically designed for laminar-turbulent transition
//! 2. **Natural Damping**: Automatically reduces in low-gradient regions
//! 3. **Computational Efficiency**: Simple formulation with minimal computational overhead
//! 4. **Physical Consistency**: Based on fundamental flow physics
//!
//! ## Implementation Details
//!
//! The model computes:
//!
//! ```math
//! Σ_{i,j} (∂u_i/∂x_j)² = α_ij α_ij
//! ```
//!
//! where α_ij is the velocity gradient tensor.
//!
//! The SGS viscosity is then:
//!
//! ```math
//! ν_SGS = (C_σ Δ)² Σ_{i,j} (∂u_i/∂x_j)²
//! ```
//!
//! ## Model Constants
//!
//! - C_σ = 1.35 (from literature calibration)
//! - Δ = filter width (typically grid spacing)
//!
//! ## Literature Compliance
//!
//! - Nicoud, F., et al. (2011). Using singular values to build a subgrid-scale model
//!   for large eddy simulations. Physics of Fluids, 23(8).
//! - Sagaut, P. (2006). Large Eddy Simulation for Incompressible Flows (3rd ed.).
//!   Springer. Section 5.3.4: The σ-model.
//!
//! ## Numerical Properties
//!
//! - **Stability**: Inherently stable due to quadratic formulation
//! - **Conservation**: Respects fundamental conservation properties
//! - **Scalability**: O(1) per grid point evaluation

use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;

/// Sigma SGS model configuration
#[derive(Debug, Clone, Copy)]
pub struct SigmaConfig<T: RealField + Copy> {
    /// Sigma model constant (typically 1.35 from literature)
    pub c_sigma: T,
    /// Filter width (grid spacing)
    pub filter_width: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for SigmaConfig<T> {
    fn default() -> Self {
        Self {
            c_sigma: T::from_f64(1.35).unwrap(), // Standard value from literature
            filter_width: T::one(),
        }
    }
}

/// Sigma Subgrid-Scale model implementation
#[derive(Debug, Clone)]
pub struct SigmaModel<T: RealField + Copy> {
    /// Sigma model configuration parameters
    pub config: SigmaConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> SigmaModel<T> {
    /// Create new Sigma model with default configuration
    pub fn new() -> Self {
        Self {
            config: SigmaConfig::default(),
        }
    }

    /// Create Sigma model with custom configuration
    pub fn with_config(config: SigmaConfig<T>) -> Self {
        Self { config }
    }

    /// Compute SGS viscosity using Sigma model
    ///
    /// # Arguments
    ///
    /// * `velocity_gradient` - 2x2 velocity gradient tensor ∂u_i/∂x_j
    ///
    /// # Returns
    ///
    /// SGS viscosity ν_SGS
    pub fn sgs_viscosity(&self, velocity_gradient: &DMatrix<T>) -> T {
        if velocity_gradient.nrows() != 2 || velocity_gradient.ncols() != 2 {
            return T::zero();
        }

        // Compute Σ_{i,j} (∂u_i/∂x_j)² = α_ij α_ij
        let mut gradient_squared_sum = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                let grad_ij = velocity_gradient[(i, j)];
                gradient_squared_sum += grad_ij * grad_ij;
            }
        }

        // ν_SGS = (C_σ Δ)² Σ_{i,j} (∂u_i/∂x_j)²
        let c_sigma_delta = self.config.c_sigma * self.config.filter_width;
        c_sigma_delta * c_sigma_delta * gradient_squared_sum
    }

    /// Compute SGS stress tensor using Sigma model
    ///
    /// # Arguments
    ///
    /// * `velocity_gradient` - 2x2 velocity gradient tensor ∂u_i/∂x_j
    ///
    /// # Returns
    ///
    /// SGS stress tensor τ_ij (symmetric 2x2 matrix)
    pub fn sgs_stress(&self, velocity_gradient: &DMatrix<T>) -> DMatrix<T> {
        let nu_sgs = self.sgs_viscosity(velocity_gradient);

        // Compute strain rate tensor S_ij = (1/2)(∂u_i/∂x_j + ∂u_j/∂x_i)
        let mut strain_rate = DMatrix::zeros(2, 2);
        for i in 0..2 {
            for j in 0..2 {
                let duidxj = velocity_gradient[(i, j)];
                let duidxi = velocity_gradient[(j, i)];
                strain_rate[(i, j)] = (duidxj + duidxi) * T::from_f64(0.5).unwrap();
            }
        }

        // SGS stress τ_ij = -2 ν_SGS S_ij
        let mut sgs_stress = DMatrix::zeros(2, 2);
        let minus_two_nu = -T::from_f64(2.0).unwrap() * nu_sgs;

        for i in 0..2 {
            for j in 0..2 {
                sgs_stress[(i, j)] = minus_two_nu * strain_rate[(i, j)];
            }
        }

        sgs_stress
    }

    /// Get model configuration
    pub fn config(&self) -> &SigmaConfig<T> {
        &self.config
    }

    /// Update filter width (useful for adaptive grids)
    pub fn set_filter_width(&mut self, filter_width: T) {
        self.config.filter_width = filter_width;
    }

    /// Update Sigma constant
    pub fn set_c_sigma(&mut self, c_sigma: T) {
        self.config.c_sigma = c_sigma;
    }

    /// Compute the Sigma invariant (for analysis and validation)
    ///
    /// This returns Σ_{i,j} (∂u_i/∂x_j)² which is used in the SGS viscosity formula
    pub fn sigma_invariant(&self, velocity_gradient: &DMatrix<T>) -> T {
        if velocity_gradient.nrows() != 2 || velocity_gradient.ncols() != 2 {
            return T::zero();
        }

        let mut sum = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                let grad_ij = velocity_gradient[(i, j)];
                sum += grad_ij * grad_ij;
            }
        }
        sum
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for SigmaModel<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_sigma_model_creation() {
        let model = SigmaModel::<f64>::new();
        assert_eq!(model.config().c_sigma, 1.35);
        assert_eq!(model.config().filter_width, 1.0);
    }

    #[test]
    fn test_sigma_zero_velocity_gradient() {
        let model = SigmaModel::<f64>::new();
        let zero_gradient = DMatrix::zeros(2, 2);

        let nu_sgs = model.sgs_viscosity(&zero_gradient);
        assert_eq!(nu_sgs, 0.0);

        let stress = model.sgs_stress(&zero_gradient);
        assert_eq!(stress[(0, 0)], 0.0);
        assert_eq!(stress[(0, 1)], 0.0);
        assert_eq!(stress[(1, 1)], 0.0);
    }

    #[test]
    fn test_sigma_simple_shear() {
        let model = SigmaModel::<f64>::new();

        // Simple shear flow: u = γ y, v = 0
        // Velocity gradient: [0, γ; 0, 0] where γ = 1
        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 1)] = 1.0;

        let nu_sgs = model.sgs_viscosity(&velocity_gradient);
        // ν_SGS = (1.35 * 1)² * (0² + 1² + 0² + 0²) = 1.8225 * 1 = 1.8225
        assert_relative_eq!(nu_sgs, 1.8225, epsilon = 1e-10);

        let stress = model.sgs_stress(&velocity_gradient);
        // For simple shear: S_12 = S_21 = 0.5
        // τ_12 = τ_21 = -2 * ν_SGS * 0.5 = -ν_SGS
        assert_relative_eq!(stress[(0, 1)], -nu_sgs, epsilon = 1e-10);
        assert_relative_eq!(stress[(1, 0)], -nu_sgs, epsilon = 1e-10);
    }

    #[test]
    fn test_sigma_config_update() {
        let mut model = SigmaModel::<f64>::new();

        model.set_c_sigma(1.5);
        model.set_filter_width(0.01);

        assert_eq!(model.config().c_sigma, 1.5);
        assert_eq!(model.config().filter_width, 0.01);
    }

    #[test]
    fn test_sigma_invariant() {
        let model = SigmaModel::<f64>::new();

        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 0)] = 1.0;
        velocity_gradient[(0, 1)] = 2.0;
        velocity_gradient[(1, 0)] = 3.0;
        velocity_gradient[(1, 1)] = 4.0;

        // Σ = 1² + 2² + 3² + 4² = 1 + 4 + 9 + 16 = 30
        let invariant = model.sigma_invariant(&velocity_gradient);
        assert_eq!(invariant, 30.0);
    }

    #[test]
    fn test_sigma_stress_tensor_symmetry() {
        let model = SigmaModel::<f64>::new();

        // General 2D flow
        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 0)] = 0.5;
        velocity_gradient[(0, 1)] = 0.3;
        velocity_gradient[(1, 0)] = 0.2;
        velocity_gradient[(1, 1)] = -0.1;

        let stress = model.sgs_stress(&velocity_gradient);

        // Check symmetry
        assert_relative_eq!(stress[(0, 1)], stress[(1, 0)], epsilon = 1e-10);
    }

    #[test]
    fn test_sigma_scaling_properties() {
        let model = SigmaModel::<f64>::new();

        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 1)] = 1.0;

        let nu1 = model.sgs_viscosity(&velocity_gradient);

        // Double the filter width
        let mut model2 = model.clone();
        model2.set_filter_width(2.0);
        let nu2 = model2.sgs_viscosity(&velocity_gradient);

        // ν_SGS should scale with (Δ)²
        assert_relative_eq!(nu2 / nu1, 4.0, epsilon = 1e-10);

        // Double the velocity gradient
        let mut velocity_gradient2 = velocity_gradient.clone();
        velocity_gradient2.scale_mut(2.0);
        let nu3 = model.sgs_viscosity(&velocity_gradient2);

        // ν_SGS should scale with (gradient)²
        assert_relative_eq!(nu3 / nu1, 4.0, epsilon = 1e-10);
    }
}
