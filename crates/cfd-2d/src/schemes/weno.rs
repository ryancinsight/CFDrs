//! Weighted Essentially Non-Oscillatory (WENO) schemes
//!
//! ## Mathematical Foundation
//!
//! WENO schemes provide high-order accuracy in smooth regions while maintaining
//! non-oscillatory behavior near discontinuities through nonlinear weighting.
//!
//! ### WENO Reconstruction
//!
//! For a 5-point stencil, WENO5 reconstructs the interface value using three
//! candidate stencils, each providing third-order accuracy:
//!
//! **Stencil 1**: {u_{j-2}, u_{j-1}, u_j} → q₁ = (2u_{j-2} - 7u_{j-1} + 11u_j)/6
//! **Stencil 2**: {u_{j-1}, u_j, u_{j+1}} → q₂ = (-u_{j-1} + 5u_j + 2u_{j+1})/6
//! **Stencil 3**: {u_j, u_{j+1}, u_{j+2}} → q₃ = (2u_j + 5u_{j+1} - u_{j+2})/6
//!
//! The final reconstruction is: u_{j+1/2} = ∑ ω_k q_k
//!
//! ## Local Truncation Error (LTE) Bounds
//!
//! ### WENO5 Scheme
//!
//! **Smooth regions**: LTE = O(Δx⁵) (fifth-order accuracy)
//!
//! **Near discontinuities**: LTE = O(Δx²) (second-order, oscillation-free)
//!
//! **LTE bound**: |τ| ≤ C Δx^p where p = 5 for smooth flows, p = 2 near shocks
//!
//! ## Stability Analysis
//!
//! ### Von Neumann Stability
//!
//! WENO schemes maintain stability similar to their underlying schemes:
//!
//! **Stability region**: CFL ≤ 1/10 for WENO5 (very restrictive)
//!
//! The nonlinear weighting provides robustness but requires small time steps.
//!
//! ### TVD Property
//!
//! WENO schemes are not strictly TVD but maintain boundedness through:
//! - Nonlinear weighting that downweights oscillatory stencils
//! - Essentially non-oscillatory behavior near discontinuities
//!
//! ## CFL Conditions
//!
//! ### WENO5: CFL ≤ 1/10 (explicit schemes)
//!
//! The restrictive CFL condition is due to:
//! - High-order accuracy requirements
//! - Nonlinear stability constraints
//! - Need to resolve small-scale features
//!
//! ## References
//!
//! - Jiang, G. S., & Shu, C. W. (1996). Efficient implementation of weighted ENO schemes.
//!   *Journal of Computational Physics*, 126(1), 202-228.
//! - Shu, C. W. (1997). Essentially non-oscillatory and weighted essentially non-oscillatory
//!   schemes for hyperbolic conservation laws. In *Advanced numerical approximation of
//!   nonlinear hyperbolic equations* (pp. 325-432). Springer.

use super::{constants, weno_constants, Grid2D, SpatialDiscretization};
use cfd_core::constants::mathematical::numeric::{THREE, TWO};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Fifth-order WENO scheme
pub struct WENO5<T: RealField + Copy> {
    epsilon: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for WENO5<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> WENO5<T> {
    /// Create new WENO5 scheme
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(constants::WENO_EPSILON).unwrap_or_else(T::zero),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Compute smoothness indicators
    fn smoothness_indicators(&self, v: &[T; 5]) -> [T; 3] {
        let coeff_13_12 =
            T::from_f64(weno_constants::WENO5_BETA_COEFF_13_12).unwrap_or_else(T::zero);
        let coeff_quarter =
            T::from_f64(weno_constants::WENO5_BETA_COEFF_QUARTER).unwrap_or_else(T::zero);
        let two = T::from_f64(TWO).unwrap_or_else(T::zero);
        let three = T::from_f64(THREE).unwrap_or_else(T::zero);
        let four = T::from_f64(weno_constants::WENO5_BETA_COEFF_FOUR).unwrap_or_else(T::zero);

        // Beta_0
        let beta0 = coeff_13_12 * (v[0] - two * v[1] + v[2]).powi(2)
            + coeff_quarter * (v[0] - four * v[1] + three * v[2]).powi(2);

        // Beta_1
        let beta1 = coeff_13_12 * (v[1] - two * v[2] + v[3]).powi(2)
            + coeff_quarter * (v[1] - v[3]).powi(2);

        // Beta_2
        let beta2 = coeff_13_12 * (v[2] - two * v[3] + v[4]).powi(2)
            + coeff_quarter * (three * v[2] - four * v[3] + v[4]).powi(2);

        [beta0, beta1, beta2]
    }

    /// Compute WENO weights
    fn weno_weights(&self, beta: &[T; 3]) -> [T; 3] {
        let d0 = T::from_f64(constants::WENO5_WEIGHTS[0]).unwrap_or_else(T::zero);
        let d1 = T::from_f64(constants::WENO5_WEIGHTS[1]).unwrap_or_else(T::zero);
        let d2 = T::from_f64(constants::WENO5_WEIGHTS[2]).unwrap_or_else(T::zero);

        let alpha0 = d0 / (self.epsilon + beta[0]).powi(2);
        let alpha1 = d1 / (self.epsilon + beta[1]).powi(2);
        let alpha2 = d2 / (self.epsilon + beta[2]).powi(2);

        let sum = alpha0 + alpha1 + alpha2;

        [alpha0 / sum, alpha1 / sum, alpha2 / sum]
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for WENO5<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // Extract stencil
        let v = [
            grid.data[(i - 2, j)],
            grid.data[(i - 1, j)],
            grid.data[(i, j)],
            grid.data[(i + 1, j)],
            grid.data[(i + 2, j)],
        ];

        // Compute smoothness indicators
        let beta = self.smoothness_indicators(&v);

        // Compute weights
        let w = self.weno_weights(&beta);

        // Compute flux
        let f0 = v[0] / T::from_f64(3.0).unwrap_or_else(T::zero)
            - T::from_f64(7.0).unwrap_or_else(T::zero) * v[1]
                / T::from_f64(6.0).unwrap_or_else(T::zero)
            + T::from_f64(11.0).unwrap_or_else(T::zero) * v[2]
                / T::from_f64(6.0).unwrap_or_else(T::zero);

        let f1 = -v[1] / T::from_f64(6.0).unwrap_or_else(T::zero)
            + T::from_f64(5.0).unwrap_or_else(T::zero) * v[2]
                / T::from_f64(6.0).unwrap_or_else(T::zero)
            + v[3] / T::from_f64(3.0).unwrap_or_else(T::zero);

        let f2 = v[2] / T::from_f64(3.0).unwrap_or_else(T::zero)
            + T::from_f64(5.0).unwrap_or_else(T::zero) * v[3]
                / T::from_f64(6.0).unwrap_or_else(T::zero)
            - v[4] / T::from_f64(6.0).unwrap_or_else(T::zero);

        (w[0] * f0 + w[1] * f1 + w[2] * f2) / grid.dx
    }

    fn order(&self) -> usize {
        5
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for WENO5 scheme: CFL ≤ 1/10
    /// WENO schemes are very dissipative and require small time steps
    fn cfl_limit(&self) -> f64 {
        0.1 // Very restrictive for stability
    }
}
