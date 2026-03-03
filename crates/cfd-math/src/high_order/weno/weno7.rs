//! WENO7 reconstruction scheme for high-order shock-capturing.
//!
//! Implements the 7th-order WENO scheme using four 4th-order stencils.

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// WENO7 reconstruction scheme for high-order shock-capturing
///
/// Implements the 7th-order WENO scheme using four 4th-order stencils.
pub struct WENO7<T: RealField + Copy + FromPrimitive> {
    /// Small parameter to avoid division by zero in smoothness indicators
    epsilon: T,
}

impl<T: RealField + Copy + FromPrimitive> WENO7<T> {
    /// Create a new WENO7 scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            epsilon: T::from_f64(1e-6).unwrap_or_else(T::zero),
        }
    }

    /// Create WENO7 with custom epsilon parameter
    #[must_use]
    pub fn with_epsilon(epsilon: T) -> Self {
        Self { epsilon }
    }

    /// Reconstruct left interface value q_{j+1/2}^- from cell averages
    ///
    /// # Arguments
    /// * `cells` - Array of cell averages [u_{j-3}, u_{j-2}, u_{j-1}, u_j, u_{j+1}, u_{j+2}, u_{j+3}]
    ///
    /// # Returns
    /// Reconstructed interface value at j+1/2
    #[must_use]
    pub fn reconstruct_left(&self, cells: &[T; 7]) -> T {
        // ENO reconstructions (4th order each)
        let q0 = self.eno4_stencil_0(cells);
        let q1 = self.eno4_stencil_1(cells);
        let q2 = self.eno4_stencil_2(cells);
        let q3 = self.eno4_stencil_3(cells);

        // Smoothness indicators
        let beta0 = self.smoothness_indicator_0(cells);
        let beta1 = self.smoothness_indicator_1(cells);
        let beta2 = self.smoothness_indicator_2(cells);
        let beta3 = self.smoothness_indicator_3(cells);

        // Linear weights for WENO7
        let c0 = T::from_f64(1.0 / 35.0).unwrap_or_else(T::one);
        let c1 = T::from_f64(12.0 / 35.0).unwrap_or_else(T::one);
        let c2 = T::from_f64(18.0 / 35.0).unwrap_or_else(T::one);
        let c3 = T::from_f64(4.0 / 35.0).unwrap_or_else(T::one);

        // Nonlinear weights
        let alpha0 = c0 / (self.epsilon + beta0).powi(2);
        let alpha1 = c1 / (self.epsilon + beta1).powi(2);
        let alpha2 = c2 / (self.epsilon + beta2).powi(2);
        let alpha3 = c3 / (self.epsilon + beta3).powi(2);

        let alpha_sum = alpha0 + alpha1 + alpha2 + alpha3;

        // Final WENO reconstruction
        let omega0 = alpha0 / alpha_sum;
        let omega1 = alpha1 / alpha_sum;
        let omega2 = alpha2 / alpha_sum;
        let omega3 = alpha3 / alpha_sum;

        omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3
    }

    /// Reconstruct right interface value q_{j+1/2}^+ from cell averages
    ///
    /// # Arguments
    /// * `cells` - Array of cell averages [u_{j-3}, u_{j-2}, u_{j-1}, u_j, u_{j+1}, u_{j+2}, u_{j+3}]
    ///
    /// # Returns
    /// Reconstructed interface value at j+1/2
    #[must_use]
    pub fn reconstruct_right(&self, cells: &[T; 7]) -> T {
        // ENO reconstructions (4th order each, right-biased)
        let q0 = self.eno4_stencil_0_right(cells);
        let q1 = self.eno4_stencil_1_right(cells);
        let q2 = self.eno4_stencil_2_right(cells);
        let q3 = self.eno4_stencil_3_right(cells);

        // Smoothness indicators (same as left reconstruction)
        let beta0 = self.smoothness_indicator_0(cells);
        let beta1 = self.smoothness_indicator_1(cells);
        let beta2 = self.smoothness_indicator_2(cells);
        let beta3 = self.smoothness_indicator_3(cells);

        // Linear weights (reversed for right reconstruction)
        let c0 = T::from_f64(4.0 / 35.0).unwrap_or_else(T::one);
        let c1 = T::from_f64(18.0 / 35.0).unwrap_or_else(T::one);
        let c2 = T::from_f64(12.0 / 35.0).unwrap_or_else(T::one);
        let c3 = T::from_f64(1.0 / 35.0).unwrap_or_else(T::one);

        // Nonlinear weights
        let alpha0 = c0 / (self.epsilon + beta0).powi(2);
        let alpha1 = c1 / (self.epsilon + beta1).powi(2);
        let alpha2 = c2 / (self.epsilon + beta2).powi(2);
        let alpha3 = c3 / (self.epsilon + beta3).powi(2);

        let alpha_sum = alpha0 + alpha1 + alpha2 + alpha3;

        // Final WENO reconstruction
        let omega0 = alpha0 / alpha_sum;
        let omega1 = alpha1 / alpha_sum;
        let omega2 = alpha2 / alpha_sum;
        let omega3 = alpha3 / alpha_sum;

        omega0 * q0 + omega1 * q1 + omega2 * q2 + omega3 * q3
    }

    /// 4th-order ENO stencil 0 (left): {u_{j-3}, u_{j-2}, u_{j-1}, u_j}
    #[must_use]
    fn eno4_stencil_0(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let neg_three = T::from_f64(-3.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let neg_twenty_three = T::from_f64(-23.0).unwrap_or_else(T::one);
        let twenty_five = T::from_f64(25.0).unwrap_or_else(T::one);

        one_twelfth
            * (neg_three * cells[0]
                + thirteen * cells[1]
                + neg_twenty_three * cells[2]
                + twenty_five * cells[3])
    }

    /// 4th-order ENO stencil 1 (left): {u_{j-2}, u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno4_stencil_1(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let neg_five = T::from_f64(-5.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let three = T::from_f64(3.0).unwrap_or_else(T::one);

        one_twelfth * (cells[1] + neg_five * cells[2] + thirteen * cells[3] + three * cells[4])
    }

    /// 4th-order ENO stencil 2 (left): {u_{j-1}, u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno4_stencil_2(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let seven = T::from_f64(7.0).unwrap_or_else(T::one);

        one_twelfth * (-cells[2] + seven * cells[3] + seven * cells[4] - cells[5])
    }

    /// 4th-order ENO stencil 3 (left): {u_j, u_{j+1}, u_{j+2}, u_{j+3}}
    #[must_use]
    fn eno4_stencil_3(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let three = T::from_f64(3.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let neg_five = T::from_f64(-5.0).unwrap_or_else(T::one);

        one_twelfth * (three * cells[3] + thirteen * cells[4] + neg_five * cells[5] + cells[6])
    }

    /// 4th-order ENO stencil 0 (right): {u_{j-3}, u_{j-2}, u_{j-1}, u_j}
    #[must_use]
    fn eno4_stencil_0_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let three = T::from_f64(3.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let neg_five = T::from_f64(-5.0).unwrap_or_else(T::one);

        one_twelfth * (cells[0] + neg_five * cells[1] + thirteen * cells[2] + three * cells[3])
    }

    /// 4th-order ENO stencil 1 (right): {u_{j-2}, u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno4_stencil_1_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let seven = T::from_f64(7.0).unwrap_or_else(T::one);

        one_twelfth * (-cells[1] + seven * cells[2] + seven * cells[3] - cells[4])
    }

    /// 4th-order ENO stencil 2 (right): {u_{j-1}, u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno4_stencil_2_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let neg_five = T::from_f64(-5.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let three = T::from_f64(3.0).unwrap_or_else(T::one);

        one_twelfth * (three * cells[2] + thirteen * cells[3] + neg_five * cells[4] + cells[5])
    }

    /// 4th-order ENO stencil 3 (right): {u_j, u_{j+1}, u_{j+2}, u_{j+3}}
    #[must_use]
    fn eno4_stencil_3_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = T::from_f64(1.0 / 12.0).unwrap_or_else(T::one);
        let neg_three = T::from_f64(-3.0).unwrap_or_else(T::one);
        let thirteen = T::from_f64(13.0).unwrap_or_else(T::one);
        let neg_twenty_three = T::from_f64(-23.0).unwrap_or_else(T::one);
        let twenty_five = T::from_f64(25.0).unwrap_or_else(T::one);

        one_twelfth
            * (twenty_five * cells[3]
                + neg_twenty_three * cells[4]
                + thirteen * cells[5]
                + neg_three * cells[6])
    }

    /// Smoothness indicator beta0 for WENO7
    #[must_use]
    fn smoothness_indicator_0(&self, cells: &[T; 7]) -> T {
        let one_240 = T::from_f64(1.0 / 240.0).unwrap_or_else(T::one);
        let c2107 = T::from_f64(2107.0).unwrap_or_else(T::one);
        let c9402 = T::from_f64(9402.0).unwrap_or_else(T::one);
        let c7042 = T::from_f64(7042.0).unwrap_or_else(T::one);
        let c1854 = T::from_f64(1854.0).unwrap_or_else(T::one);
        let c11003 = T::from_f64(11003.0).unwrap_or_else(T::one);
        let c17246 = T::from_f64(17246.0).unwrap_or_else(T::one);
        let c4642 = T::from_f64(4642.0).unwrap_or_else(T::one);
        let c7043 = T::from_f64(7043.0).unwrap_or_else(T::one);
        let c3882 = T::from_f64(3882.0).unwrap_or_else(T::one);
        let c547 = T::from_f64(547.0).unwrap_or_else(T::one);

        let u0 = cells[0];
        let u1 = cells[1];
        let u2 = cells[2];
        let u3 = cells[3];

        one_240
            * (u0 * (c2107 * u0 - c9402 * u1 + c7042 * u2 - c1854 * u3)
                + u1 * (c11003 * u1 - c17246 * u2 + c4642 * u3)
                + u2 * (c7043 * u2 - c3882 * u3)
                + c547 * u3 * u3)
    }

    /// Smoothness indicator beta1 for WENO7
    #[must_use]
    fn smoothness_indicator_1(&self, cells: &[T; 7]) -> T {
        let one_240 = T::from_f64(1.0 / 240.0).unwrap_or_else(T::one);
        let c547 = T::from_f64(547.0).unwrap_or_else(T::one);
        let c2522 = T::from_f64(2522.0).unwrap_or_else(T::one);
        let c1922 = T::from_f64(1922.0).unwrap_or_else(T::one);
        let c494 = T::from_f64(494.0).unwrap_or_else(T::one);
        let c3443 = T::from_f64(3443.0).unwrap_or_else(T::one);
        let c5966 = T::from_f64(5966.0).unwrap_or_else(T::one);
        let c1602 = T::from_f64(1602.0).unwrap_or_else(T::one);
        let c2843 = T::from_f64(2843.0).unwrap_or_else(T::one);
        let c1642 = T::from_f64(1642.0).unwrap_or_else(T::one);
        let c267 = T::from_f64(267.0).unwrap_or_else(T::one);

        let u1 = cells[1];
        let u2 = cells[2];
        let u3 = cells[3];
        let u4 = cells[4];

        one_240
            * (u1 * (c547 * u1 - c2522 * u2 + c1922 * u3 - c494 * u4)
                + u2 * (c3443 * u2 - c5966 * u3 + c1602 * u4)
                + u3 * (c2843 * u3 - c1642 * u4)
                + c267 * u4 * u4)
    }

    /// Smoothness indicator beta2 for WENO7
    #[must_use]
    fn smoothness_indicator_2(&self, cells: &[T; 7]) -> T {
        let one_240 = T::from_f64(1.0 / 240.0).unwrap_or_else(T::one);
        let c267 = T::from_f64(267.0).unwrap_or_else(T::one);
        let c1642 = T::from_f64(1642.0).unwrap_or_else(T::one);
        let c1602 = T::from_f64(1602.0).unwrap_or_else(T::one);
        let c494 = T::from_f64(494.0).unwrap_or_else(T::one);
        let c2843 = T::from_f64(2843.0).unwrap_or_else(T::one);
        let c5966 = T::from_f64(5966.0).unwrap_or_else(T::one);
        let c1922 = T::from_f64(1922.0).unwrap_or_else(T::one);
        let c3443 = T::from_f64(3443.0).unwrap_or_else(T::one);
        let c2522 = T::from_f64(2522.0).unwrap_or_else(T::one);
        let c547 = T::from_f64(547.0).unwrap_or_else(T::one);

        let u2 = cells[2];
        let u3 = cells[3];
        let u4 = cells[4];
        let u5 = cells[5];

        one_240
            * (u2 * (c267 * u2 - c1642 * u3 + c1602 * u4 - c494 * u5)
                + u3 * (c2843 * u3 - c5966 * u4 + c1922 * u5)
                + u4 * (c3443 * u4 - c2522 * u5)
                + c547 * u5 * u5)
    }

    /// Smoothness indicator beta3 for WENO7
    #[must_use]
    fn smoothness_indicator_3(&self, cells: &[T; 7]) -> T {
        let one_240 = T::from_f64(1.0 / 240.0).unwrap_or_else(T::one);
        let c547 = T::from_f64(547.0).unwrap_or_else(T::one);
        let c3882 = T::from_f64(3882.0).unwrap_or_else(T::one);
        let c4642 = T::from_f64(4642.0).unwrap_or_else(T::one);
        let c1854 = T::from_f64(1854.0).unwrap_or_else(T::one);
        let c7043 = T::from_f64(7043.0).unwrap_or_else(T::one);
        let c17246 = T::from_f64(17246.0).unwrap_or_else(T::one);
        let c7042 = T::from_f64(7042.0).unwrap_or_else(T::one);
        let c11003 = T::from_f64(11003.0).unwrap_or_else(T::one);
        let c9402 = T::from_f64(9402.0).unwrap_or_else(T::one);
        let c2107 = T::from_f64(2107.0).unwrap_or_else(T::one);

        let u3 = cells[3];
        let u4 = cells[4];
        let u5 = cells[5];
        let u6 = cells[6];

        one_240
            * (u3 * (c547 * u3 - c3882 * u4 + c4642 * u5 - c1854 * u6)
                + u4 * (c7043 * u4 - c17246 * u5 + c7042 * u6)
                + u5 * (c11003 * u5 - c9402 * u6)
                + c2107 * u6 * u6)
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for WENO7<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_weno7_smooth_function() {
        let weno = WENO7::<f64>::new();
        let dx = 0.1;
        let x_j = 1.0;
        let x_interface = x_j + 0.5 * dx;

        let avg = |x: f64| ((x - 0.5 * dx).cos() - (x + 0.5 * dx).cos()) / dx;

        let cells = [
            avg(x_j - 3.0 * dx),
            avg(x_j - 2.0 * dx),
            avg(x_j - dx),
            avg(x_j),
            avg(x_j + dx),
            avg(x_j + 2.0 * dx),
            avg(x_j + 3.0 * dx),
        ];

        let reconstructed = weno.reconstruct_left(&cells);
        let analytical = x_interface.sin();

        assert!((reconstructed - analytical).abs() < 1e-5);
    }

    #[test]
    fn test_weno7_smoothness_indicators() {
        let weno = WENO7::<f64>::new();

        let constant_cells = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        assert_relative_eq!(weno.smoothness_indicator_0(&constant_cells), 0.0);
        assert_relative_eq!(weno.smoothness_indicator_1(&constant_cells), 0.0);
        assert_relative_eq!(weno.smoothness_indicator_2(&constant_cells), 0.0);
        assert_relative_eq!(weno.smoothness_indicator_3(&constant_cells), 0.0);
    }
}
