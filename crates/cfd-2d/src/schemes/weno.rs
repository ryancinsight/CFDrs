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
use cfd_core::physics::constants::mathematical::numeric::{THREE, TWO};
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

/// Ninth-order WENO scheme
///
/// ## Mathematical Foundation
///
/// WENO9 provides ninth-order accuracy in smooth regions using an 11-point stencil.
/// The scheme uses 5 candidate stencils, each providing fifth-order accuracy.
///
/// ### Stencil Configuration
/// WENO9 uses 5 candidate stencils over 11 points:
/// - Stencil 0: {u_{j-5}, ..., u_{j}}     (5th-order ENO)
/// - Stencil 1: {u_{j-4}, ..., u_{j+1}}   (5th-order ENO)
/// - Stencil 2: {u_{j-3}, ..., u_{j+2}}   (5th-order ENO)
/// - Stencil 3: {u_{j-2}, ..., u_{j+3}}   (5th-order ENO)
/// - Stencil 4: {u_{j-1}, ..., u_{j+4}}   (5th-order ENO)
///
/// ### Local Truncation Error
/// - **Smooth regions**: LTE = O(Δx⁹) (ninth-order accuracy)
/// - **Near discontinuities**: LTE = O(Δx⁵) (fifth-order, oscillation-free)
///
/// ### CFL Condition
/// WENO9 requires CFL ≤ 1/18 due to extreme high-order accuracy requirements
/// and nonlinear stability constraints.
pub struct WENO9<T: RealField + Copy> {
    epsilon: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> Default for WENO9<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> WENO9<T> {
    /// Create new WENO9 scheme
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(constants::WENO_EPSILON).unwrap_or_else(T::zero),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Compute smoothness indicators for WENO9 (5 indicators for 11-point stencil)
    fn smoothness_indicators(&self, v: &[T; 11]) -> [T; 5] {
        // WENO9 smoothness indicators based on Jiang-Shu formulation
        // These are the optimized coefficients for 9th-order accuracy

        let mut beta = [T::zero(); 5];

        // Beta_0 (stencil 0: u[j-5..j])
        beta[0] = T::from_f64(0.0015308084989341916).unwrap_or_else(T::zero)
            * (v[0] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[1]
                + T::from_f64(5.0).unwrap_or_else(T::zero) * v[2])
                .powi(2)
            + T::from_f64(0.002_740_988_421_902_865).unwrap_or_else(T::zero)
                * (v[0] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[1]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[2]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    + v[4])
                    .powi(2)
            + T::from_f64(0.031254897785245544).unwrap_or_else(T::zero)
                * (v[0] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[1]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[2]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    + T::from_f64(5.0).unwrap_or_else(T::zero) * v[5])
                    .powi(2);

        // Beta_1 (stencil 1: u[j-4..j+1])
        beta[1] = T::from_f64(0.0015308084989341916).unwrap_or_else(T::zero)
            * (v[1] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[2]
                + T::from_f64(5.0).unwrap_or_else(T::zero) * v[3])
                .powi(2)
            + T::from_f64(0.002_740_988_421_902_865).unwrap_or_else(T::zero)
                * (v[1] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[2]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    + v[5])
                    .powi(2)
            + T::from_f64(0.031254897785245544).unwrap_or_else(T::zero)
                * (v[1] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[2]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    + T::from_f64(5.0).unwrap_or_else(T::zero) * v[6])
                    .powi(2);

        // Beta_2 (stencil 2: u[j-3..j+2])
        beta[2] = T::from_f64(0.0015308084989341916).unwrap_or_else(T::zero)
            * (v[2] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                + T::from_f64(5.0).unwrap_or_else(T::zero) * v[4])
                .powi(2)
            + T::from_f64(0.002_740_988_421_902_865).unwrap_or_else(T::zero)
                * (v[2] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    + v[6])
                    .powi(2)
            + T::from_f64(0.031254897785245544).unwrap_or_else(T::zero)
                * (v[2] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[3]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[6]
                    + T::from_f64(5.0).unwrap_or_else(T::zero) * v[7])
                    .powi(2);

        // Beta_3 (stencil 3: u[j-2..j+3])
        beta[3] = T::from_f64(0.0015308084989341916).unwrap_or_else(T::zero)
            * (v[3] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                + T::from_f64(5.0).unwrap_or_else(T::zero) * v[5])
                .powi(2)
            + T::from_f64(0.002_740_988_421_902_865).unwrap_or_else(T::zero)
                * (v[3] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[6]
                    + v[7])
                    .powi(2)
            + T::from_f64(0.031254897785245544).unwrap_or_else(T::zero)
                * (v[3] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[4]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[6]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[7]
                    + T::from_f64(5.0).unwrap_or_else(T::zero) * v[8])
                    .powi(2);

        // Beta_4 (stencil 4: u[j-1..j+4])
        beta[4] = T::from_f64(0.0015308084989341916).unwrap_or_else(T::zero)
            * (v[4] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                + T::from_f64(5.0).unwrap_or_else(T::zero) * v[6])
                .powi(2)
            + T::from_f64(0.002_740_988_421_902_865).unwrap_or_else(T::zero)
                * (v[4] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[6]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[7]
                    + v[8])
                    .powi(2)
            + T::from_f64(0.031254897785245544).unwrap_or_else(T::zero)
                * (v[4] - T::from_f64(4.0).unwrap_or_else(T::zero) * v[5]
                    + T::from_f64(4.0).unwrap_or_else(T::zero) * v[6]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[7]
                    - T::from_f64(4.0).unwrap_or_else(T::zero) * v[8]
                    + T::from_f64(5.0).unwrap_or_else(T::zero) * v[9])
                    .powi(2);

        beta
    }

    /// Compute WENO9 weights using optimized coefficients
    fn weno_weights(&self, beta: &[T; 5]) -> [T; 5] {
        // Optimized weights for WENO9 (Henrick et al. 2005)
        let d = [
            T::from_f64(weno_constants::WENO9_LINEAR_WEIGHTS[0]).unwrap_or_else(T::zero),
            T::from_f64(weno_constants::WENO9_LINEAR_WEIGHTS[1]).unwrap_or_else(T::zero),
            T::from_f64(weno_constants::WENO9_LINEAR_WEIGHTS[2]).unwrap_or_else(T::zero),
            T::from_f64(weno_constants::WENO9_LINEAR_WEIGHTS[3]).unwrap_or_else(T::zero),
            T::from_f64(weno_constants::WENO9_LINEAR_WEIGHTS[4]).unwrap_or_else(T::zero),
        ];

        let mut alpha = [T::zero(); 5];
        for i in 0..5 {
            alpha[i] = d[i] / (self.epsilon + beta[i]).powi(2);
        }

        let sum: T = alpha.iter().copied().sum();
        let mut weights = [T::zero(); 5];
        for i in 0..5 {
            weights[i] = alpha[i] / sum;
        }

        weights
    }
}

impl<T: RealField + Copy + FromPrimitive + std::iter::Sum> SpatialDiscretization<T> for WENO9<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // Extract 11-point stencil (requires boundary checking in real implementation)
        let v = [
            grid.data[(i - 5, j)],
            grid.data[(i - 4, j)],
            grid.data[(i - 3, j)],
            grid.data[(i - 2, j)],
            grid.data[(i - 1, j)],
            grid.data[(i, j)],
            grid.data[(i + 1, j)],
            grid.data[(i + 2, j)],
            grid.data[(i + 3, j)],
            grid.data[(i + 4, j)],
            grid.data[(i + 5, j)],
        ];

        // Compute smoothness indicators
        let beta = self.smoothness_indicators(&v);

        // Compute weights
        let w = self.weno_weights(&beta);

        // Compute reconstructed flux using 5 candidate stencils
        let mut flux = T::zero();
        let denom =
            T::from_f64(weno_constants::WENO9_STENCIL_DENOM).unwrap_or_else(T::zero);

        for k in 0..5 {
            let mut q_k = T::zero();
            // Stencil k uses points from v[1+k] to v[1+k+4]
            // v indices correspond to: v[5] is u_i
            // k=0: v[1]..v[5] (u_{i-4}..u_i)
            // ...
            // k=4: v[5]..v[9] (u_i..u_{i+4})
            for j in 0..5 {
                let coeff =
                    T::from_f64(weno_constants::WENO9_STENCIL_COEFFS[k][j])
                        .unwrap_or_else(T::zero);
                q_k += coeff * v[1 + k + j];
            }
            q_k /= denom;
            flux += w[k] * q_k;
        }

        flux / grid.dx
    }

    fn order(&self) -> usize {
        9
    }

    fn is_conservative(&self) -> bool {
        true
    }

    /// CFL limit for WENO9 scheme: CFL ≤ 1/18
    /// Extremely restrictive due to 9th-order accuracy requirements
    fn cfl_limit(&self) -> f64 {
        1.0 / 18.0 // Extremely restrictive for stability
    }
}
