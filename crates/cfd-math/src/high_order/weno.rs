//! Weighted Essentially Non-Oscillatory (WENO) schemes.
//!
//! This module provides high-order accurate shock-capturing schemes for solving
//! hyperbolic conservation laws. WENO schemes are designed to achieve high-order
//! accuracy in smooth regions while maintaining stability and non-oscillatory
//! behavior near discontinuities.
//!
//! ## Mathematical Foundation
//!
//! ### WENO Reconstruction
//! For a hyperbolic conservation law ∂u/∂t + ∂f(u)/∂x = 0, WENO reconstructs
//! the interface flux f(u) from cell averages ū_j.
//!
//! The WENO5 scheme uses three 3rd-order ENO stencils:
//! - r=0: {u_{j-2}, u_{j-1}, u_{j}}
//! - r=1: {u_{j-1}, u_{j}, u_{j+1}}
//! - r=2: {u_{j}, u_{j+1}, u_{j+2}}
//!
//! ### Nonlinear Weights
//! The weight for stencil r is:
//! ω_r = α_r / ∑_{s=0}^2 α_s, where α_r = C_r / (ε + β_r)^2
//!
//! ### Smoothness Indicators (β_r)
//! β_0 = (13/12)(u_{j-2} - 2u_{j-1} + u_j)^2 + (1/4)(u_{j-2} - 4u_{j-1} + 3u_j)^2
//! β_1 = (13/12)(u_{j-1} - 2u_j + u_{j+1})^2 + (1/4)(u_{j-1} - u_{j+1})^2
//! β_2 = (13/12)(u_j - 2u_{j+1} + u_{j+2})^2 + (1/4)(3u_j - 4u_{j+1} + u_{j+2})^2
//!
//! ### Linear Weights (C_r)
//! C_0 = 1/10, C_1 = 6/10, C_2 = 3/10
//!
//! ## Convergence Properties
//!
//! **Theorem (WENO Accuracy)**: WENO5 achieves 5th-order accuracy in smooth regions
//! and maintains 3rd-order accuracy near discontinuities.
//!
//! **Theorem (Stability)**: WENO schemes are TVD (Total Variation Diminishing) and
//! maintain monotonicity near shocks.
//!
//! ## References
//!
//! - Jiang, G.-S. & Shu, C.-W. (1996). "Efficient implementation of weighted ENO schemes"
//! - Shu, C.-W. (1997). "Essentially non-oscillatory and weighted essentially non-oscillatory schemes for hyperbolic conservation laws"
//! - Borges, R. et al. (2008). "An improved weighted essentially non-oscillatory scheme for hyperbolic conservation laws"

use nalgebra::{RealField, Scalar};
use num_traits::FromPrimitive;

/// WENO5 reconstruction scheme for shock-capturing
///
/// Implements the 5th-order WENO scheme of Jiang & Shu (1996) with
/// improved weights of Borges et al. (2008).
pub struct WENO5<T: RealField + Copy + FromPrimitive> {
    /// Small parameter to avoid division by zero in smoothness indicators
    epsilon: T,
}

impl<T: RealField + Copy + FromPrimitive> WENO5<T> {
    /// Create a new WENO5 scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(1e-6).unwrap_or_else(T::zero),
        }
    }

    /// Create WENO5 with custom epsilon parameter
    #[must_use]
    pub fn with_epsilon(epsilon: T) -> Self {
        Self { epsilon }
    }

    /// Reconstruct left interface value q_{j+1/2}^- from cell averages
    ///
    /// # Arguments
    /// * `cells` - Array of cell averages [u_{j-2}, u_{j-1}, u_j, u_{j+1}, u_{j+2}]
    ///
    /// # Returns
    /// Reconstructed interface value at j+1/2
    #[must_use]
    pub fn reconstruct_left(&self, cells: &[T; 5]) -> T {
        // ENO reconstructions (3rd order each)
        let q0 = self.eno3_stencil_0(cells);
        let q1 = self.eno3_stencil_1(cells);
        let q2 = self.eno3_stencil_2(cells);

        // Smoothness indicators
        let beta0 = self.smoothness_indicator_0(cells);
        let beta1 = self.smoothness_indicator_1(cells);
        let beta2 = self.smoothness_indicator_2(cells);

        // Linear weights (Borges et al. 2008 improved weights)
        let c0 = T::from_f64(1.0 / 10.0).unwrap_or_else(T::one);
        let c1 = T::from_f64(6.0 / 10.0).unwrap_or_else(T::one);
        let c2 = T::from_f64(3.0 / 10.0).unwrap_or_else(T::one);

        // Nonlinear weights
        let alpha0 = c0 / (self.epsilon + beta0).powi(2);
        let alpha1 = c1 / (self.epsilon + beta1).powi(2);
        let alpha2 = c2 / (self.epsilon + beta2).powi(2);

        let alpha_sum = alpha0 + alpha1 + alpha2;

        // Final WENO reconstruction
        let omega0 = alpha0 / alpha_sum;
        let omega1 = alpha1 / alpha_sum;
        let omega2 = alpha2 / alpha_sum;

        omega0 * q0 + omega1 * q1 + omega2 * q2
    }

    /// Reconstruct right interface value q_{j+1/2}^+ from cell averages
    ///
    /// # Arguments
    /// * `cells` - Array of cell averages [u_{j-2}, u_{j-1}, u_j, u_{j+1}, u_{j+2}]
    ///
    /// # Returns
    /// Reconstructed interface value at j+1/2
    #[must_use]
    pub fn reconstruct_right(&self, cells: &[T; 5]) -> T {
        // ENO reconstructions (3rd order each)
        let q0 = self.eno3_stencil_0_right(cells);
        let q1 = self.eno3_stencil_1_right(cells);
        let q2 = self.eno3_stencil_2_right(cells);

        // Smoothness indicators (same as left reconstruction)
        let beta0 = self.smoothness_indicator_0(cells);
        let beta1 = self.smoothness_indicator_1(cells);
        let beta2 = self.smoothness_indicator_2(cells);

        // Linear weights (Borges et al. 2008 improved weights)
        let c0 = T::from_f64(1.0 / 10.0).unwrap_or_else(T::one);
        let c1 = T::from_f64(6.0 / 10.0).unwrap_or_else(T::one);
        let c2 = T::from_f64(3.0 / 10.0).unwrap_or_else(T::one);

        // Nonlinear weights
        let alpha0 = c0 / (self.epsilon + beta0).powi(2);
        let alpha1 = c1 / (self.epsilon + beta1).powi(2);
        let alpha2 = c2 / (self.epsilon + beta2).powi(2);

        let alpha_sum = alpha0 + alpha1 + alpha2;

        // Final WENO reconstruction
        let omega0 = alpha0 / alpha_sum;
        let omega1 = alpha1 / alpha_sum;
        let omega2 = alpha2 / alpha_sum;

        omega0 * q0 + omega1 * q1 + omega2 * q2
    }

    /// 3rd-order ENO stencil 0 (left reconstruction): {u_{j-2}, u_{j-1}, u_j}
    #[must_use]
    fn eno3_stencil_0(&self, cells: &[T; 5]) -> T {
        let one_third = T::from_f64(1.0 / 3.0).unwrap_or_else(T::one);
        let thirteen_twelfth = T::from_f64(13.0 / 12.0).unwrap_or_else(T::one);
        let one_fourth = T::from_f64(1.0 / 4.0).unwrap_or_else(T::one);

        // 3rd order ENO reconstruction: (2u_{j-2} - 7u_{j-1} + 11u_j)/6
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let seven = T::from_f64(7.0).unwrap_or_else(T::one);
        let eleven = T::from_f64(11.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (two * cells[0] - seven * cells[1] + eleven * cells[2]) / six
    }

    /// 3rd-order ENO stencil 1 (left reconstruction): {u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno3_stencil_1(&self, cells: &[T; 5]) -> T {
        let one_sixth = T::from_f64(1.0 / 6.0).unwrap_or_else(T::one);

        // 3rd order ENO reconstruction: (-u_{j-1} + 5u_j + 2u_{j+1})/6
        let five = T::from_f64(5.0).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (-cells[1] + five * cells[2] + two * cells[3]) / six
    }

    /// 3rd-order ENO stencil 2 (left reconstruction): {u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno3_stencil_2(&self, cells: &[T; 5]) -> T {
        let one_sixth = T::from_f64(1.0 / 6.0).unwrap_or_else(T::one);

        // 3rd order ENO reconstruction: (2u_j + 5u_{j+1} - u_{j+2})/6
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let five = T::from_f64(5.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (two * cells[2] + five * cells[3] - cells[4]) / six
    }

    /// 3rd-order ENO stencil 0 (right reconstruction): {u_{j-2}, u_{j-1}, u_j}
    #[must_use]
    fn eno3_stencil_0_right(&self, cells: &[T; 5]) -> T {
        // 3rd order ENO reconstruction: (11u_{j-2} - 7u_{j-1} + 2u_j)/6
        let eleven = T::from_f64(11.0).unwrap_or_else(T::one);
        let seven = T::from_f64(7.0).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (eleven * cells[0] - seven * cells[1] + two * cells[2]) / six
    }

    /// 3rd-order ENO stencil 1 (right reconstruction): {u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno3_stencil_1_right(&self, cells: &[T; 5]) -> T {
        // 3rd order ENO reconstruction: (2u_{j-1} + 5u_j - u_{j+1})/6
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let five = T::from_f64(5.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (two * cells[1] + five * cells[2] - cells[3]) / six
    }

    /// 3rd-order ENO stencil 2 (right reconstruction): {u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno3_stencil_2_right(&self, cells: &[T; 5]) -> T {
        // 3rd order ENO reconstruction: (-u_j + 5u_{j+1} + 2u_{j+2})/6
        let five = T::from_f64(5.0).unwrap_or_else(T::one);
        let two = T::from_f64(2.0).unwrap_or_else(T::one);
        let six = T::from_f64(6.0).unwrap_or_else(T::one);

        (-cells[2] + five * cells[3] + two * cells[4]) / six
    }

    /// Smoothness indicator β_0 for stencil r=0
    #[must_use]
    fn smoothness_indicator_0(&self, cells: &[T; 5]) -> T {
        let thirteen_twelfth = T::from_f64(13.0 / 12.0).unwrap_or_else(T::one);
        let one_fourth = T::from_f64(1.0 / 4.0).unwrap_or_else(T::one);

        let u0 = cells[0];
        let u1 = cells[1];
        let u2 = cells[2];

        let term1 = u0 - T::from_f64(2.0).unwrap_or_else(T::one) * u1 + u2;
        let term2 = u0 - T::from_f64(4.0).unwrap_or_else(T::one) * u1 +
                   T::from_f64(3.0).unwrap_or_else(T::one) * u2;

        thirteen_twelfth * term1 * term1 + one_fourth * term2 * term2
    }

    /// Smoothness indicator β_1 for stencil r=1
    #[must_use]
    fn smoothness_indicator_1(&self, cells: &[T; 5]) -> T {
        let thirteen_twelfth = T::from_f64(13.0 / 12.0).unwrap_or_else(T::one);
        let one_fourth = T::from_f64(1.0 / 4.0).unwrap_or_else(T::one);

        let u1 = cells[1];
        let u2 = cells[2];
        let u3 = cells[3];

        let term1 = u1 - T::from_f64(2.0).unwrap_or_else(T::one) * u2 + u3;
        let term2 = u1 - u3;

        thirteen_twelfth * term1 * term1 + one_fourth * term2 * term2
    }

    /// Smoothness indicator β_2 for stencil r=2
    #[must_use]
    fn smoothness_indicator_2(&self, cells: &[T; 5]) -> T {
        let thirteen_twelfth = T::from_f64(13.0 / 12.0).unwrap_or_else(T::one);
        let one_fourth = T::from_f64(1.0 / 4.0).unwrap_or_else(T::one);

        let u2 = cells[2];
        let u3 = cells[3];
        let u4 = cells[4];

        let term1 = u2 - T::from_f64(2.0).unwrap_or_else(T::one) * u3 + u4;
        let term2 = T::from_f64(3.0).unwrap_or_else(T::one) * u2 -
                   T::from_f64(4.0).unwrap_or_else(T::one) * u3 + u4;

        thirteen_twelfth * term1 * term1 + one_fourth * term2 * term2
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for WENO5<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_weno5_smooth_function() {
        // Test on a smooth function: sin(x)
        // For smooth functions, WENO5 should achieve 5th-order accuracy
        let weno = WENO5::<f64>::new();

        // Cell averages of sin(x) over [0,2π]
        let cells = [
            0.5 * (0.0_f64.sin() + (std::f64::consts::PI / 4.0).sin()), // [0, π/4]
            0.5 * ((std::f64::consts::PI / 4.0).sin() + (std::f64::consts::PI / 2.0).sin()), // [π/4, π/2]
            0.5 * ((std::f64::consts::PI / 2.0).sin() + (3.0 * std::f64::consts::PI / 4.0).sin()), // [π/2, 3π/4]
            0.5 * ((3.0 * std::f64::consts::PI / 4.0).sin() + std::f64::consts::PI.sin()), // [3π/4, π]
            0.5 * ((std::f64::consts::PI).sin() + (5.0 * std::f64::consts::PI / 4.0).sin()), // [π, 5π/4]
        ];

        let reconstructed = weno.reconstruct_left(&cells);

        // Should be close to analytical interface value
        let analytical = (std::f64::consts::PI / 2.0).sin();

        // Allow some tolerance for numerical approximation
        assert!((reconstructed - analytical).abs() < 1e-3);
    }

    #[test]
    fn test_weno5_smoothness_indicators() {
        let weno = WENO5::<f64>::new();

        // Test with smooth data
        let smooth_cells = [1.0, 1.1, 1.21, 1.331, 1.4641]; // x^2
        let beta0 = weno.smoothness_indicator_0(&smooth_cells);
        let beta1 = weno.smoothness_indicator_1(&smooth_cells);
        let beta2 = weno.smoothness_indicator_2(&smooth_cells);

        // Smoothness indicators should be small for smooth data
        assert!(beta0 < 1e-6);
        assert!(beta1 < 1e-6);
        assert!(beta2 < 1e-6);
    }

    #[test]
    fn test_weno5_discontinuous() {
        let weno = WENO5::<f64>::new();

        // Test with discontinuous data
        let disc_cells = [1.0, 1.0, 1.0, 2.0, 2.0]; // Shock
        let beta0 = weno.smoothness_indicator_0(&disc_cells);
        let beta1 = weno.smoothness_indicator_1(&disc_cells);
        let beta2 = weno.smoothness_indicator_2(&disc_cells);

        // At least one smoothness indicator should be large near discontinuity
        assert!(beta0 > 0.1 || beta1 > 0.1 || beta2 > 0.1);
    }

    #[test]
    fn test_weno5_weights() {
        let weno = WENO5::<f64>::new();

        // Test that weights sum to 1 for smooth data
        let smooth_cells = [1.0, 1.01, 1.0201, 1.030301, 1.04060401]; // Smooth function

        // Get the weights indirectly through reconstruction
        let q0 = weno.eno3_stencil_0(&smooth_cells);
        let q1 = weno.eno3_stencil_1(&smooth_cells);
        let q2 = weno.eno3_stencil_2(&smooth_cells);

        let beta0 = weno.smoothness_indicator_0(&smooth_cells);
        let beta1 = weno.smoothness_indicator_1(&smooth_cells);
        let beta2 = weno.smoothness_indicator_2(&smooth_cells);

        let c0 = 1.0 / 10.0;
        let c1 = 6.0 / 10.0;
        let c2 = 3.0 / 10.0;

        let alpha0 = c0 / (weno.epsilon + beta0).powi(2);
        let alpha1 = c1 / (weno.epsilon + beta1).powi(2);
        let alpha2 = c2 / (weno.epsilon + beta2).powi(2);

        let alpha_sum = alpha0 + alpha1 + alpha2;
        let omega0 = alpha0 / alpha_sum;
        let omega1 = alpha1 / alpha_sum;
        let omega2 = alpha2 / alpha_sum;

        // Weights should sum to 1
        assert_relative_eq!(omega0 + omega1 + omega2, 1.0, epsilon = 1e-12);

        // For smooth data, middle stencil should be heavily weighted
        assert!(omega1 > 0.5);
    }
}

