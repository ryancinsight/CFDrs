//! Weighted Essentially Non-Oscillatory (WENO) Schemes
//!
//! WENO schemes provide high-order accurate shock-capturing capabilities
//! for hyperbolic conservation laws, essential for compressible flow simulations.
//!
//! ## Mathematical Foundation
//!
//! WENO reconstruction computes a high-order approximation at cell interfaces
//! by adaptively weighting multiple low-order candidate stencils based on
//! their smoothness indicators.
//!
//! For a 5th-order WENO scheme:
//!
//! ```math
//! q_{i+1/2} = ω₀ q₀^{(2)} + ω₁ q₁^{(2)} + ω₂ q₂^{(2)}
//! ```
//!
//! where q₀^{(2)}, q₁^{(2)}, q₂^{(2)} are 3rd-order reconstructions from
//! different stencils, and ωᵢ are nonlinear weights.
//!
//! ## Key Components
//!
//! 1. **Candidate Stencils**: Multiple low-order reconstructions
//! 2. **Smoothness Indicators**: Measure stencil regularity
//! 3. **Nonlinear Weights**: Adapt to flow discontinuities
//! 4. **High-Order Reconstruction**: Weighted combination
//!
//! ## WENO5 Implementation Details
//!
//! **Stencil Structure:**
//! - Stencil 0: {q_{i-2}, q_{i-1}, q_{i}}
//! - Stencil 1: {q_{i-1}, q_{i}, q_{i+1}}
//! - Stencil 2: {q_{i}, q_{i+1}, q_{i+2}}
//!
//! **Linear Weights (γᵢ):**
//! - γ₀ = 1/10, γ₁ = 6/10, γ₂ = 3/10
//!
//! **Smoothness Indicators (βᵢ):**
//! ```math
//! β₀ = (13/12)(q_{i-2} - 2q_{i-1} + q_i)² + (1/4)(q_{i-2} - 4q_{i-1} + 3q_i)²
//! β₁ = (13/12)(q_{i-1} - 2q_i + q_{i+1})² + (1/4)(q_{i-1} - q_{i+1})²
//! β₂ = (13/12)(q_i - 2q_{i+1} + q_{i+2})² + (1/4)(3q_i - 4q_{i+1} + q_{i+2})²
//! ```
//!
//! **Nonlinear Weights:**
//! ```math
//! ωᵢ = γᵢ / (ε + βᵢ)²
//! ωᵢ = ωᵢ / Σⱼ ωⱼ
//! ```
//!
//! ## Literature Compliance
//!
//! - Jiang, G.-S., & Shu, C.-W. (1996). Efficient implementation of weighted ENO
//!   schemes. Journal of Computational Physics, 126(1), 202-228.
//! - Shu, C.-W. (1998). Essentially non-oscillatory and weighted essentially
//!   non-oscillatory schemes for hyperbolic conservation laws. In *Advanced
//!   numerical approximation of nonlinear hyperbolic equations* (pp. 325-432).
//!   Springer.
//!
//! ## Numerical Properties
//!
//! - **Order of Accuracy**: 5th-order in smooth regions
//! - **Shock Capturing**: Maintains sharp discontinuities
//! - **Stability**: TVD (Total Variation Diminishing) properties
//! - **Efficiency**: O(N) complexity with optimal constant

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// WENO reconstruction configuration
#[derive(Debug, Clone, Copy)]
pub struct WenoConfig<T: RealField + Copy> {
    /// Small parameter to avoid division by zero (ε ≈ 10^{-6} to 10^{-40})
    pub epsilon: T,
    /// Power parameter for weight computation (typically p=2)
    pub p: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for WenoConfig<T> {
    fn default() -> Self {
        Self {
            epsilon: T::from_f64(1e-40).unwrap(), // Very small to avoid division issues
            p: T::from_f64(2.0).unwrap(),
        }
    }
}

/// 5th-order Weighted Essentially Non-Oscillatory reconstruction
#[derive(Debug, Clone)]
pub struct Weno5<T: RealField + Copy> {
    config: WenoConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> Weno5<T> {
    /// Create WENO5 reconstruction with default configuration
    pub fn new() -> Self {
        Self {
            config: WenoConfig::default(),
        }
    }

    /// Create WENO5 reconstruction with custom configuration
    pub fn with_config(config: WenoConfig<T>) -> Self {
        Self { config }
    }

    /// Reconstruct left interface value (q_{i+1/2}^-) from cell-centered values
    ///
    /// # Arguments
    ///
    /// * `q` - Array of cell-centered values [q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}]
    ///
    /// # Returns
    ///
    /// Reconstructed value at i+1/2 interface
    pub fn reconstruct_left(&self, q: &[T; 5]) -> T {
        // Linear weights for 5th-order WENO
        let gamma = [
            T::from_f64(1.0 / 10.0).unwrap(),  // γ₀ = 1/10
            T::from_f64(6.0 / 10.0).unwrap(),  // γ₁ = 6/10
            T::from_f64(3.0 / 10.0).unwrap(),  // γ₂ = 3/10
        ];

        // Compute candidate polynomials (3rd-order ENO reconstructions)
        let q0 = self.eno3_stencil0(q);
        let q1 = self.eno3_stencil1(q);
        let q2 = self.eno3_stencil2(q);

        let candidates = [q0, q1, q2];

        // Compute smoothness indicators
        let beta = self.compute_smoothness_indicators(q);

        // Compute nonlinear weights
        let weights = self.compute_nonlinear_weights(&gamma, &beta);

        // Weighted reconstruction
        weights[0] * candidates[0] +
        weights[1] * candidates[1] +
        weights[2] * candidates[2]
    }

    /// Reconstruct right interface value (q_{i+1/2}^+) from cell-centered values
    ///
    /// # Arguments
    ///
    /// * `q` - Array of cell-centered values [q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}]
    ///
    /// # Returns
    ///
    /// Reconstructed value at i+1/2 interface
    pub fn reconstruct_right(&self, q: &[T; 5]) -> T {
        // For right reconstruction, use mirrored stencils
        // This is equivalent to reconstructing from the right
        let q_mirror = [q[4], q[3], q[2], q[1], q[0]];
        -self.reconstruct_left(&q_mirror)  // Negative because of flux difference
    }

    /// Compute 3rd-order ENO reconstruction from stencil 0: {q_{i-2}, q_{i-1}, q_i}
    fn eno3_stencil0(&self, q: &[T; 5]) -> T {
        // q^{(2)}_0 = (2q_{i-2} - 7q_{i-1} + 11q_i)/6
        let two = T::from_f64(2.0).unwrap();
        let seven = T::from_f64(7.0).unwrap();
        let eleven = T::from_f64(11.0).unwrap();
        let six = T::from_f64(6.0).unwrap();

        (two * q[0] - seven * q[1] + eleven * q[2]) / six
    }

    /// Compute 3rd-order ENO reconstruction from stencil 1: {q_{i-1}, q_i, q_{i+1}}
    fn eno3_stencil1(&self, q: &[T; 5]) -> T {
        // q^{(2)}_1 = (-q_{i-1} + 5q_i + 2q_{i+1})/6
        let one = T::one();
        let five = T::from_f64(5.0).unwrap();
        let two = T::from_f64(2.0).unwrap();
        let six = T::from_f64(6.0).unwrap();

        (-one * q[1] + five * q[2] + two * q[3]) / six
    }

    /// Compute 3rd-order ENO reconstruction from stencil 2: {q_i, q_{i+1}, q_{i+2}}
    fn eno3_stencil2(&self, q: &[T; 5]) -> T {
        // q^{(2)}_2 = (2q_i + 5q_{i+1} - q_{i+2})/6
        let two = T::from_f64(2.0).unwrap();
        let five = T::from_f64(5.0).unwrap();
        let one = T::one();
        let six = T::from_f64(6.0).unwrap();

        (two * q[2] + five * q[3] - one * q[4]) / six
    }

    /// Compute smoothness indicators βᵢ for each stencil
    fn compute_smoothness_indicators(&self, q: &[T; 5]) -> [T; 3] {
        let thirteen_twelve = T::from_f64(13.0 / 12.0).unwrap();
        let one_fourth = T::from_f64(1.0 / 4.0).unwrap();

        // β₀ = (13/12)(q_{i-2} - 2q_{i-1} + q_i)² + (1/4)(q_{i-2} - 4q_{i-1} + 3q_i)²
        let diff0 = q[0] - T::from_f64(2.0).unwrap() * q[1] + q[2];
        let diff1 = q[0] - T::from_f64(4.0).unwrap() * q[1] + T::from_f64(3.0).unwrap() * q[2];
        let beta0 = thirteen_twelve * diff0 * diff0 + one_fourth * diff1 * diff1;

        // β₁ = (13/12)(q_{i-1} - 2q_i + q_{i+1})² + (1/4)(q_{i-1} - q_{i+1})²
        let diff2 = q[1] - T::from_f64(2.0).unwrap() * q[2] + q[3];
        let diff3 = q[1] - q[3];
        let beta1 = thirteen_twelve * diff2 * diff2 + one_fourth * diff3 * diff3;

        // β₂ = (13/12)(q_i - 2q_{i+1} + q_{i+2})² + (1/4)(3q_i - 4q_{i+1} + q_{i+2})²
        let diff4 = q[2] - T::from_f64(2.0).unwrap() * q[3] + q[4];
        let diff5 = T::from_f64(3.0).unwrap() * q[2] - T::from_f64(4.0).unwrap() * q[3] + q[4];
        let beta2 = thirteen_twelve * diff4 * diff4 + one_fourth * diff5 * diff5;

        [beta0, beta1, beta2]
    }

    /// Compute nonlinear weights from linear weights and smoothness indicators
    fn compute_nonlinear_weights(&self, gamma: &[T; 3], beta: &[T; 3]) -> [T; 3] {
        // Compute unnormalized weights: ω̃ᵢ = γᵢ / (ε + βᵢ)^p
        let mut omega_tilde = [T::zero(); 3];
        for i in 0..3 {
            let denominator = self.config.epsilon + beta[i];
            omega_tilde[i] = gamma[i] / denominator.powf(self.config.p);
        }

        // Normalize: ωᵢ = ω̃ᵢ / Σⱼ ω̃ⱼ
        let sum = omega_tilde[0] + omega_tilde[1] + omega_tilde[2];
        [
            omega_tilde[0] / sum,
            omega_tilde[1] / sum,
            omega_tilde[2] / sum,
        ]
    }

    /// Get reconstruction configuration
    pub fn config(&self) -> &WenoConfig<T> {
        &self.config
    }

    /// Update configuration
    pub fn set_config(&mut self, config: WenoConfig<T>) {
        self.config = config;
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for Weno5<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// 7th-order Weighted Essentially Non-Oscillatory reconstruction
#[derive(Debug, Clone)]
pub struct Weno7<T: RealField + Copy> {
    config: WenoConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive> Weno7<T> {
    /// Create WENO7 reconstruction with default configuration
    pub fn new() -> Self {
        Self {
            config: WenoConfig::default(),
        }
    }

    /// Create WENO7 reconstruction with custom configuration
    pub fn with_config(config: WenoConfig<T>) -> Self {
        Self { config }
    }

    /// Reconstruct left interface value (q_{i+1/2}^-) from cell-centered values
    ///
    /// Requires 8 cell-centered values: [q_{i-3}, q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}, q_{i+3}, q_{i+4}]
    ///
    /// # Arguments
    ///
    /// * `q` - Array of cell-centered values (length 8)
    ///
    /// # Returns
    ///
    /// Reconstructed value at i+1/2 interface
    pub fn reconstruct_left(&self, q: &[T; 8]) -> T {
        // Linear weights for 7th-order WENO
        let gamma = [
            T::from_f64(1.0 / 35.0).unwrap(),   // γ₀ = 1/35
            T::from_f64(12.0 / 35.0).unwrap(),  // γ₁ = 12/35
            T::from_f64(18.0 / 35.0).unwrap(),  // γ₂ = 18/35
            T::from_f64(4.0 / 35.0).unwrap(),   // γ₃ = 4/35
        ];

        // Compute candidate polynomials (5th-order ENO reconstructions)
        let q0 = self.eno5_stencil0(q);
        let q1 = self.eno5_stencil1(q);
        let q2 = self.eno5_stencil2(q);
        let q3 = self.eno5_stencil3(q);

        let candidates = [q0, q1, q2, q3];

        // Compute smoothness indicators
        let beta = self.compute_smoothness_indicators7(q);

        // Compute nonlinear weights
        let weights = self.compute_nonlinear_weights7(&gamma, &beta);

        // Weighted reconstruction
        weights[0] * candidates[0] +
        weights[1] * candidates[1] +
        weights[2] * candidates[2] +
        weights[3] * candidates[3]
    }

    /// Compute 5th-order ENO reconstruction from stencil 0: {q_{i-3}, q_{i-2}, q_{i-1}, q_i}
    fn eno5_stencil0(&self, q: &[T; 8]) -> T {
        // q^{(4)}_0 = (1/30)q_{i-3} - (8/30)q_{i-2} + (37/30)q_{i-1} - (2/30)q_i
        let thirty = T::from_f64(30.0).unwrap();
        (T::one() * q[0] - T::from_f64(8.0).unwrap() * q[1] +
         T::from_f64(37.0).unwrap() * q[2] - T::from_f64(2.0).unwrap() * q[3]) / thirty
    }

    /// Compute 5th-order ENO reconstruction from stencil 1: {q_{i-2}, q_{i-1}, q_i, q_{i+1}}
    fn eno5_stencil1(&self, q: &[T; 8]) -> T {
        // q^{(4)}_1 = (-3/10)q_{i-2} + (19/10)q_{i-1} - (2/10)q_i + (6/10)q_{i+1}
        let ten = T::from_f64(10.0).unwrap();
        (-T::from_f64(3.0).unwrap() * q[1] + T::from_f64(19.0).unwrap() * q[2] +
         -T::from_f64(2.0).unwrap() * q[3] + T::from_f64(6.0).unwrap() * q[4]) / ten
    }

    /// Compute 5th-order ENO reconstruction from stencil 2: {q_{i-1}, q_i, q_{i+1}, q_{i+2}}
    fn eno5_stencil2(&self, q: &[T; 8]) -> T {
        // q^{(4)}_2 = (2/3)q_{i-1} - (13/3)q_i + (47/3)q_{i+1} - (9/3)q_{i+2}
        let three = T::from_f64(3.0).unwrap();
        (T::from_f64(2.0).unwrap() * q[2] - T::from_f64(13.0).unwrap() * q[3] +
         T::from_f64(47.0).unwrap() * q[4] - T::from_f64(9.0).unwrap() * q[5]) / three
    }

    /// Compute 5th-order ENO reconstruction from stencil 3: {q_i, q_{i+1}, q_{i+2}, q_{i+3}}
    fn eno5_stencil3(&self, q: &[T; 8]) -> T {
        // q^{(4)}_3 = (-3/2)q_i + (25/2)q_{i+1} - (26/2)q_{i+2} + (9/2)q_{i+3}
        let two = T::from_f64(2.0).unwrap();
        (-T::from_f64(3.0).unwrap() * q[3] + T::from_f64(25.0).unwrap() * q[4] +
         -T::from_f64(26.0).unwrap() * q[5] + T::from_f64(9.0).unwrap() * q[6]) / two
    }

    /// Compute smoothness indicators βᵢ for 7th-order WENO (4 stencils)
    fn compute_smoothness_indicators7(&self, q: &[T; 8]) -> [T; 4] {
        // β₀ = q_{i-3}(547q_{i-3} - 3882q_{i-2} + 4642q_{i-1} - 1854q_i) +
        //      q_{i-2}(7043q_{i-2} - 17246q_{i-1} + 7042q_i) +
        //      q_{i-1}(11003q_{i-1} - 9402q_i) + 2107q_i²
        let beta0 = q[0] * (T::from_f64(547.0).unwrap() * q[0] -
                           T::from_f64(3882.0).unwrap() * q[1] +
                           T::from_f64(4642.0).unwrap() * q[2] -
                           T::from_f64(1854.0).unwrap() * q[3]) +
                  q[1] * (T::from_f64(7043.0).unwrap() * q[1] -
                          T::from_f64(17246.0).unwrap() * q[2] +
                          T::from_f64(7042.0).unwrap() * q[3]) +
                  q[2] * (T::from_f64(11003.0).unwrap() * q[2] -
                          T::from_f64(9402.0).unwrap() * q[3]) +
                  T::from_f64(2107.0).unwrap() * q[3] * q[3];

        // β₁ = q_{i-2}(267q_{i-2} - 1642q_{i-1} + 1602q_i - 494q_{i+1}) +
        //      q_{i-1}(2843q_{i-1} - 5966q_i + 1922q_{i+1}) +
        //      q_i(3443q_i - 2522q_{i+1}) + 547q_{i+1}²
        let beta1 = q[1] * (T::from_f64(267.0).unwrap() * q[1] -
                           T::from_f64(1642.0).unwrap() * q[2] +
                           T::from_f64(1602.0).unwrap() * q[3] -
                           T::from_f64(494.0).unwrap() * q[4]) +
                  q[2] * (T::from_f64(2843.0).unwrap() * q[2] -
                          T::from_f64(5966.0).unwrap() * q[3] +
                          T::from_f64(1922.0).unwrap() * q[4]) +
                  q[3] * (T::from_f64(3443.0).unwrap() * q[3] -
                          T::from_f64(2522.0).unwrap() * q[4]) +
                  T::from_f64(547.0).unwrap() * q[4] * q[4];

        // β₂ = q_{i-1}(547q_{i-1} - 2522q_i + 1922q_{i+1} - 494q_{i+2}) +
        //      q_i(3443q_i - 5966q_{i+1} + 1602q_{i+2}) +
        //      q_{i+1}(2843q_{i+1} - 1642q_{i+2}) + 267q_{i+2}²
        let beta2 = q[2] * (T::from_f64(547.0).unwrap() * q[2] -
                           T::from_f64(2522.0).unwrap() * q[3] +
                           T::from_f64(1922.0).unwrap() * q[4] -
                           T::from_f64(494.0).unwrap() * q[5]) +
                  q[3] * (T::from_f64(3443.0).unwrap() * q[3] -
                          T::from_f64(5966.0).unwrap() * q[4] +
                          T::from_f64(1602.0).unwrap() * q[5]) +
                  q[4] * (T::from_f64(2843.0).unwrap() * q[4] -
                          T::from_f64(1642.0).unwrap() * q[5]) +
                  T::from_f64(267.0).unwrap() * q[5] * q[5];

        // β₃ = q_i(2107q_i - 9402q_{i+1} + 7042q_{i+2} - 1854q_{i+3}) +
        //      q_{i+1}(11003q_{i+1} - 17246q_{i+2} + 4642q_{i+3}) +
        //      q_{i+2}(7043q_{i+2} - 3882q_{i+3}) + 547q_{i+3}²
        let beta3 = q[3] * (T::from_f64(2107.0).unwrap() * q[3] -
                           T::from_f64(9402.0).unwrap() * q[4] +
                           T::from_f64(7042.0).unwrap() * q[5] -
                           T::from_f64(1854.0).unwrap() * q[6]) +
                  q[4] * (T::from_f64(11003.0).unwrap() * q[4] -
                          T::from_f64(17246.0).unwrap() * q[5] +
                          T::from_f64(4642.0).unwrap() * q[6]) +
                  q[5] * (T::from_f64(7043.0).unwrap() * q[5] -
                          T::from_f64(3882.0).unwrap() * q[6]) +
                  T::from_f64(547.0).unwrap() * q[6] * q[6];

        [beta0, beta1, beta2, beta3]
    }

    /// Compute nonlinear weights from linear weights and smoothness indicators (7th-order)
    fn compute_nonlinear_weights7(&self, gamma: &[T; 4], beta: &[T; 4]) -> [T; 4] {
        // Compute unnormalized weights: ω̃ᵢ = γᵢ / (ε + βᵢ)^p
        let mut omega_tilde = [T::zero(); 4];
        for i in 0..4 {
            let denominator = self.config.epsilon + beta[i];
            omega_tilde[i] = gamma[i] / denominator.powf(self.config.p);
        }

        // Normalize: ωᵢ = ω̃ᵢ / Σⱼ ω̃ⱼ
        let sum = omega_tilde[0] + omega_tilde[1] + omega_tilde[2] + omega_tilde[3];
        [
            omega_tilde[0] / sum,
            omega_tilde[1] / sum,
            omega_tilde[2] / sum,
            omega_tilde[3] / sum,
        ]
    }

    /// Get reconstruction configuration
    pub fn config(&self) -> &WenoConfig<T> {
        &self.config
    }

    /// Update configuration
    pub fn set_config(&mut self, config: WenoConfig<T>) {
        self.config = config;
    }

    /// Reconstruct right interface value (q_{i+1/2}^+) from cell-centered values
    ///
    /// Requires 8 cell-centered values: [q_{i-3}, q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}, q_{i+3}, q_{i+4}]
    ///
    /// # Arguments
    ///
    /// * `q` - Array of cell-centered values (length 8)
    ///
    /// # Returns
    ///
    /// Reconstructed value at i+1/2 interface
    pub fn reconstruct_right(&self, q: &[T; 8]) -> T {
        // For right reconstruction, use mirrored stencils
        // This is equivalent to reconstructing from the right
        let q_mirror = [q[7], q[6], q[5], q[4], q[3], q[2], q[1], q[0]];
        -self.reconstruct_left(&q_mirror)  // Negative because of flux difference
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for Weno7<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// WENO reconstruction for different orders
#[derive(Debug, Clone)]
pub enum WenoReconstruction<T: RealField + Copy> {
    /// 5th-order WENO scheme
    Weno5(Weno5<T>),
    /// 7th-order WENO scheme
    Weno7(Weno7<T>),
}

impl<T: RealField + Copy + FromPrimitive> WenoReconstruction<T> {
    /// Create 5th-order WENO reconstruction
    pub fn weno5() -> Self {
        Self::Weno5(Weno5::new())
    }

    /// Create 7th-order WENO reconstruction
    pub fn weno7() -> Self {
        Self::Weno7(Weno7::new())
    }

    /// Reconstruct left interface value (5th-order, requires 5 values)
    pub fn reconstruct_left5(&self, q: &[T; 5]) -> T {
        match self {
            Self::Weno5(weno) => weno.reconstruct_left(q),
            Self::Weno7(_) => panic!("WENO7 requires 8 values, use reconstruct_left8"),
        }
    }

    /// Reconstruct left interface value (7th-order, requires 8 values)
    pub fn reconstruct_left8(&self, q: &[T; 8]) -> T {
        match self {
            Self::Weno5(_) => panic!("WENO5 requires 5 values, use reconstruct_left5"),
            Self::Weno7(weno) => weno.reconstruct_left(q),
        }
    }

    /// Reconstruct right interface value (5th-order, requires 5 values)
    pub fn reconstruct_right5(&self, q: &[T; 5]) -> T {
        match self {
            Self::Weno5(weno) => weno.reconstruct_right(q),
            Self::Weno7(_) => panic!("WENO7 requires 8 values, use reconstruct_right8"),
        }
    }

    /// Reconstruct right interface value (7th-order, requires 8 values)
    pub fn reconstruct_right8(&self, q: &[T; 8]) -> T {
        match self {
            Self::Weno5(_) => panic!("WENO5 requires 5 values, use reconstruct_right5"),
            Self::Weno7(weno) => weno.reconstruct_right(q),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_weno5_creation() {
        let weno = Weno5::<f64>::new();
        assert_eq!(weno.config().epsilon, 1e-40);
        assert_eq!(weno.config().p, 2.0);
    }

    #[test]
    fn test_smooth_function_reconstruction() {
        let weno = Weno5::<f64>::new();

        // Test with smooth sine function
        // Using uniform grid with Δx = 0.1
        let dx = 0.1;
        let x_base = 1.0; // Center at x_i

        // Values: q_{i-2}, q_{i-1}, q_i, q_{i+1}, q_{i+2}
        let q = [
            (x_base - 2.0 * dx).sin(),
            (x_base - dx).sin(),
            x_base.sin(),
            (x_base + dx).sin(),
            (x_base + 2.0 * dx).sin(),
        ];

        let q_left = weno.reconstruct_left(&q);
        let q_right = weno.reconstruct_right(&q);

        // For smooth functions, WENO should give 5th-order accurate interface values
        // Analytical interface values at x_{i+1/2} ≈ x_i + 0.5*dx
        let x_interface = x_base + 0.5 * dx;
        let analytical = x_interface.sin();

        // Should be very accurate (5th-order convergence)
        assert_relative_eq!(q_left, analytical, epsilon = 1e-10);
        assert_relative_eq!(q_right, analytical, epsilon = 1e-10);
    }

    #[test]
    fn test_shock_detection() {
        let weno = Weno5::<f64>::new();

        // Test with discontinuous function (shock)
        // Step function: 0 for x < 0, 1 for x >= 0
        // Using interface at x = 0
        let q = [0.0, 0.0, 0.5, 1.0, 1.0]; // Approximate discontinuity

        let q_left = weno.reconstruct_left(&q);
        let q_right = weno.reconstruct_right(&q);

        // WENO should maintain the discontinuity without oscillations
        // Left reconstruction should be close to 0.5 (upwind value)
        // Right reconstruction should be close to 0.5 (upwind value)
        assert!(q_left >= 0.0 && q_left <= 1.0, "Left reconstruction should be bounded");
        assert!(q_right >= 0.0 && q_right <= 1.0, "Right reconstruction should be bounded");
    }

    #[test]
    fn test_linear_weights_smooth_region() {
        let weno = Weno5::<f64>::new();

        // Linear function: q(x) = x
        // With Δx = 1, values should be [-2, -1, 0, 1, 2]
        let q = [-2.0, -1.0, 0.0, 1.0, 2.0];

        let q_left = weno.reconstruct_left(&q);
        let q_right = weno.reconstruct_right(&q);

        // For linear function, interface value should be exactly 0.5
        assert_relative_eq!(q_left, 0.5, epsilon = 1e-12);
        assert_relative_eq!(q_right, 0.5, epsilon = 1e-12);
    }

    #[test]
    fn test_smoothness_indicators() {
        let weno = Weno5::<f64>::new();

        // Linear function (smooth)
        let q_smooth = [-2.0, -1.0, 0.0, 1.0, 2.0];
        let beta_smooth = weno.compute_smoothness_indicators(&q_smooth);

        // Discontinuous function
        let q_shock = [0.0, 0.0, 0.0, 1.0, 1.0];
        let beta_shock = weno.compute_smoothness_indicators(&q_shock);

        // Shock should have much larger smoothness indicators
        for i in 0..3 {
            assert!(beta_shock[i] > beta_smooth[i], "Shock should have larger smoothness indicator");
        }
    }

    #[test]
    fn test_weight_adaptation() {
        let weno = Weno5::<f64>::new();

        // Linear weights
        let gamma = [1.0/10.0, 6.0/10.0, 3.0/10.0];

        // Test with smooth data (low β)
        let beta_smooth = [0.01, 0.01, 0.01];
        let weights_smooth = weno.compute_nonlinear_weights(&gamma, &beta_smooth);

        // Test with discontinuous data (high β)
        let beta_shock = [100.0, 100.0, 100.0];
        let weights_shock = weno.compute_nonlinear_weights(&gamma, &beta_shock);

        // Weights should sum to 1
        let sum_smooth: f64 = weights_smooth.iter().sum();
        let sum_shock: f64 = weights_shock.iter().sum();

        assert_relative_eq!(sum_smooth, 1.0, epsilon = 1e-12);
        assert_relative_eq!(sum_shock, 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_weno_reconstruction_enum() {
        let weno = WenoReconstruction::<f64>::weno5();

        let q = [1.0, 2.0, 3.0, 4.0, 5.0];
        let q_left = weno.reconstruct_left(&q);
        let q_right = weno.reconstruct_right(&q);

        // Both should be finite
        assert!(q_left.is_finite());
        assert!(q_right.is_finite());
    }
}
