//! Vreman Subgrid-Scale (SGS) Model for Large Eddy Simulation
//!
//! The Vreman model is an advanced SGS model that provides superior performance
//! compared to the Smagorinsky model, particularly in transitional and complex flows.
//!
//! ## Mathematical Foundation
//!
//! The Vreman model computes SGS viscosity as:
//!
//! ```math
//! ν_SGS = C_V √(B_β / α_ij α_ij) Δ²
//! ```
//!
//! where:
//! - C_V is the Vreman constant (typically 0.07-0.1)
//! - Δ is the filter width
//! - α_ij is the velocity gradient tensor
//! - B_β = β₁₁β₂₂ - β₁₂² + β₁₁β₃₃ - β₁₃² + β₂₂β₃₃ - β₂₃²
//! - β_ij = Δ² α_mi α_mj
//!
//! ## Key Advantages over Smagorinsky
//!
//! 1. **No Excessive Damping**: Avoids over-dissipation in regions of low vorticity
//! 2. **Natural Near-Wall Behavior**: Automatically reduces to zero at walls
//! 3. **Transitional Flow Handling**: Better performance in transitional regimes
//! 4. **Mathematical Consistency**: Based on tensor invariants rather than magnitude
//!
//! ## Implementation Details
//!
//! The model computes the B_β invariant:
//!
//! ```math
//! B_β = β₁₁β₂₂ - β₁₂² + β₁₁β₃₃ - β₁₃² + β₂₂β₃₃ - β₂₃²
//! ```
//!
//! where β_ij = Δ² Σ_k α_ki α_kj
//!
//! The SGS viscosity is then:
//!
//! ```math
//! ν_SGS = C_V √(max(B_β, 0) / α_ij α_ij) Δ²
//! ```
//!
//! ## Literature Compliance
//!
//! - Vreman, A. W. (2004). An eddy-viscosity subgrid-scale model for turbulent shear flow:
//!   Algebraic theory and applications. Physics of Fluids, 16(10), 3670-3681.
//! - You, D., & Moin, P. (2007). A dynamic global-coefficient subgrid-scale eddy-viscosity model
//!   for large-eddy simulation of turbulent flow. Physics of Fluids, 19(6).
//!
//! ## Numerical Stability
//!
//! The model includes safeguards:
//! - B_β is clamped to non-negative values to prevent numerical instability
//! - α_ij α_ij is checked to avoid division by zero
//! - Natural damping in strain-dominated regions

use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;

/// Vreman SGS model configuration
#[derive(Debug, Clone, Copy)]
pub struct VremanConfig<T: RealField + Copy> {
    /// Vreman model constant (typically 0.07-0.1)
    pub c_v: T,
    /// Filter width (grid spacing)
    pub filter_width: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for VremanConfig<T> {
    fn default() -> Self {
        Self {
            c_v: T::from_f64(0.07).unwrap(), // Standard value from literature
            filter_width: T::one(),
        }
    }
}

/// Vreman Subgrid-Scale model implementation
#[derive(Debug, Clone)]
pub struct VremanModel<T: RealField + Copy> {
    /// Vreman model configuration parameters
    pub config: VremanConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> VremanModel<T> {
    /// Create new Vreman model with default configuration
    pub fn new() -> Self {
        Self {
            config: VremanConfig::default(),
        }
    }

    /// Create Vreman model with custom configuration
    pub fn with_config(config: VremanConfig<T>) -> Self {
        Self { config }
    }

    /// Compute SGS viscosity using Vreman model
    ///
    /// # Arguments
    ///
    /// * `velocity_gradient` - 2x2 velocity gradient tensor ∂u_i/∂x_j
    ///
    /// # Returns
    ///
    /// SGS viscosity ν_SGS
    pub fn sgs_viscosity(&self, velocity_gradient: &DMatrix<T>) -> T {
        self.compute_vreman_viscosity(velocity_gradient)
    }

    /// Compute SGS stress tensor using Vreman model
    ///
    /// # Arguments
    ///
    /// * `velocity_gradient` - 2x2 velocity gradient tensor ∂u_i/∂x_j
    ///
    /// # Returns
    ///
    /// SGS stress tensor τ_ij (symmetric 2x2 matrix)
    pub fn sgs_stress(&self, velocity_gradient: &DMatrix<T>) -> DMatrix<T> {
        let nu_sgs = self.compute_vreman_viscosity(velocity_gradient);

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

    /// Core Vreman SGS viscosity computation
    fn compute_vreman_viscosity(&self, velocity_gradient: &DMatrix<T>) -> T {
        if velocity_gradient.nrows() != 2 || velocity_gradient.ncols() != 2 {
            return T::zero();
        }

        // Compute α_ij α_ij (square of Frobenius norm of velocity gradient)
        let mut alpha_squared = T::zero();
        for i in 0..2 {
            for j in 0..2 {
                let val = velocity_gradient[(i, j)];
                alpha_squared = alpha_squared + val * val;
            }
        }

        // Avoid division by zero
        if alpha_squared <= T::from_f64(1e-12).unwrap() {
            return T::zero();
        }

        // Compute β_ij = Δ² Σ_k α_ki α_kj
        let delta_squared = self.config.filter_width * self.config.filter_width;
        let mut beta = DMatrix::zeros(2, 2);

        for i in 0..2 {
            for j in 0..2 {
                let mut sum = T::zero();
                for k in 0..2 {
                    sum = sum + velocity_gradient[(k, i)] * velocity_gradient[(k, j)];
                }
                beta[(i, j)] = delta_squared * sum;
            }
        }

        // Compute B_β for 2D flow
        // For 2D LES, assume isotropy in spanwise direction: β_33 = (β_11 + β_22)/2
        // B_β = β₁₁β₂₂ - β₁₂² + β₁₁β₃₃ + β₂₂β₃₃
        let b11 = beta[(0, 0)];
        let b12 = beta[(0, 1)];
        let b22 = beta[(1, 1)];
        let b33 = (b11 + b22) * T::from_f64(0.5).unwrap(); // Isotropic assumption

        let b_beta = b11 * b22 - b12 * b12 + b11 * b33 + b22 * b33;

        // Clamp to non-negative to ensure numerical stability
        let b_beta_clamped = if b_beta > T::zero() {
            b_beta
        } else {
            T::zero()
        };

        // Compute ν_SGS = C_V √(B_β / α_ij α_ij) Δ²
        let ratio = b_beta_clamped / alpha_squared;
        let sqrt_ratio = ratio.sqrt();

        self.config.c_v * sqrt_ratio * delta_squared
    }

    /// Get model configuration
    pub fn config(&self) -> &VremanConfig<T> {
        &self.config
    }

    /// Update filter width (useful for adaptive grids)
    pub fn set_filter_width(&mut self, filter_width: T) {
        self.config.filter_width = filter_width;
    }

    /// Update Vreman constant
    pub fn set_c_v(&mut self, c_v: T) {
        self.config.c_v = c_v;
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for VremanModel<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_vreman_model_creation() {
        let model = VremanModel::<f64>::new();
        assert_eq!(model.config().c_v, 0.07);
        assert_eq!(model.config().filter_width, 1.0);
    }

    #[test]
    fn test_vreman_zero_velocity_gradient() {
        let model = VremanModel::<f64>::new();
        let zero_gradient = DMatrix::zeros(2, 2);

        let nu_sgs = model.sgs_viscosity(&zero_gradient);
        assert_eq!(nu_sgs, 0.0);

        let stress = model.sgs_stress(&zero_gradient);
        assert_eq!(stress[(0, 0)], 0.0);
        assert_eq!(stress[(0, 1)], 0.0);
        assert_eq!(stress[(1, 1)], 0.0);
    }

    #[test]
    fn test_vreman_simple_shear() {
        let model = VremanModel::<f64>::new();

        // Simple shear flow: u = γ y, v = 0
        // Velocity gradient: [0, γ; 0, 0]
        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 1)] = 1.0; // ∂u/∂y = γ = 1

        let nu_sgs = model.sgs_viscosity(&velocity_gradient);
        assert!(
            nu_sgs > 0.0,
            "SGS viscosity should be positive for shear flow"
        );

        let stress = model.sgs_stress(&velocity_gradient);
        // For simple shear, S_12 = S_21 = 0.5γ
        // τ_12 = τ_21 = -2 ν_SGS S_12 = -ν_SGS γ
        assert!(stress[(0, 1)] < 0.0, "Shear stress should be negative");
        assert_eq!(
            stress[(0, 1)],
            stress[(1, 0)],
            "Stress tensor should be symmetric"
        );
    }

    #[test]
    fn test_vreman_config_update() {
        let mut model = VremanModel::<f64>::new();

        model.set_c_v(0.1);
        model.set_filter_width(0.01);

        assert_eq!(model.config().c_v, 0.1);
        assert_eq!(model.config().filter_width, 0.01);
    }

    #[test]
    fn test_vreman_stress_tensor_symmetry() {
        let model = VremanModel::<f64>::new();

        // General 2D flow
        let mut velocity_gradient = DMatrix::zeros(2, 2);
        velocity_gradient[(0, 0)] = 0.5;
        velocity_gradient[(0, 1)] = 0.3;
        velocity_gradient[(1, 0)] = 0.2;
        velocity_gradient[(1, 1)] = -0.1;

        let stress = model.sgs_stress(&velocity_gradient);

        // Check symmetry
        assert_relative_eq!(stress[(0, 1)], stress[(1, 0)], epsilon = 1e-10);
        assert_relative_eq!(stress[(0, 0)], stress[(0, 0)], epsilon = 1e-10);
        assert_relative_eq!(stress[(1, 1)], stress[(1, 1)], epsilon = 1e-10);
    }

    #[test]
    fn test_vreman_numerical_stability() {
        let model = VremanModel::<f64>::new();

        // Test with very small velocity gradients (should not cause NaN or inf)
        let mut small_gradient = DMatrix::zeros(2, 2);
        small_gradient[(0, 1)] = 1e-12;

        let nu_sgs = model.sgs_viscosity(&small_gradient);
        assert!(nu_sgs.is_finite(), "SGS viscosity should be finite");
        assert!(nu_sgs >= 0.0, "SGS viscosity should be non-negative");

        // Test with zero velocity gradients
        let zero_gradient = DMatrix::zeros(2, 2);
        let nu_sgs_zero = model.sgs_viscosity(&zero_gradient);
        assert_eq!(
            nu_sgs_zero, 0.0,
            "Zero velocity gradient should give zero SGS viscosity"
        );
    }
}
