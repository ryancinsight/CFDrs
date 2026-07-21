//! WENO7 reconstruction scheme for high-order shock-capturing.
//!
//! Implements the 7th-order WENO scheme using four 4th-order stencils.

use eunomia::{FloatElement, RealField};

use super::{from_f64, squared};

/// WENO7 reconstruction scheme for high-order shock-capturing
///
/// Implements the 7th-order WENO scheme using four 4th-order stencils.
pub struct WENO7<T: RealField + Copy + FloatElement> {
    /// Small parameter to avoid division by zero in smoothness indicators
    epsilon: T,
}

impl<T: RealField + Copy + FloatElement> WENO7<T> {
    /// Create a new WENO7 scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            epsilon: from_f64::<T>(1e-6),
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
        let c0 = from_f64::<T>(1.0 / 35.0);
        let c1 = from_f64::<T>(12.0 / 35.0);
        let c2 = from_f64::<T>(18.0 / 35.0);
        let c3 = from_f64::<T>(4.0 / 35.0);

        // Nonlinear weights
        let alpha0 = c0 / squared(self.epsilon + beta0);
        let alpha1 = c1 / squared(self.epsilon + beta1);
        let alpha2 = c2 / squared(self.epsilon + beta2);
        let alpha3 = c3 / squared(self.epsilon + beta3);

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
        let c0 = from_f64::<T>(4.0 / 35.0);
        let c1 = from_f64::<T>(18.0 / 35.0);
        let c2 = from_f64::<T>(12.0 / 35.0);
        let c3 = from_f64::<T>(1.0 / 35.0);

        // Nonlinear weights
        let alpha0 = c0 / squared(self.epsilon + beta0);
        let alpha1 = c1 / squared(self.epsilon + beta1);
        let alpha2 = c2 / squared(self.epsilon + beta2);
        let alpha3 = c3 / squared(self.epsilon + beta3);

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
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let neg_three = from_f64::<T>(-3.0);
        let thirteen = from_f64::<T>(13.0);
        let neg_twenty_three = from_f64::<T>(-23.0);
        let twenty_five = from_f64::<T>(25.0);

        one_twelfth
            * (neg_three * cells[0]
                + thirteen * cells[1]
                + neg_twenty_three * cells[2]
                + twenty_five * cells[3])
    }

    /// 4th-order ENO stencil 1 (left): {u_{j-2}, u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno4_stencil_1(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let neg_five = from_f64::<T>(-5.0);
        let thirteen = from_f64::<T>(13.0);
        let three = from_f64::<T>(3.0);

        one_twelfth * (cells[1] + neg_five * cells[2] + thirteen * cells[3] + three * cells[4])
    }

    /// 4th-order ENO stencil 2 (left): {u_{j-1}, u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno4_stencil_2(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let seven = from_f64::<T>(7.0);

        one_twelfth * (-cells[2] + seven * cells[3] + seven * cells[4] - cells[5])
    }

    /// 4th-order ENO stencil 3 (left): {u_j, u_{j+1}, u_{j+2}, u_{j+3}}
    #[must_use]
    fn eno4_stencil_3(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let three = from_f64::<T>(3.0);
        let thirteen = from_f64::<T>(13.0);
        let neg_five = from_f64::<T>(-5.0);

        one_twelfth * (three * cells[3] + thirteen * cells[4] + neg_five * cells[5] + cells[6])
    }

    /// 4th-order ENO stencil 0 (right): {u_{j-3}, u_{j-2}, u_{j-1}, u_j}
    #[must_use]
    fn eno4_stencil_0_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let three = from_f64::<T>(3.0);
        let thirteen = from_f64::<T>(13.0);
        let neg_five = from_f64::<T>(-5.0);

        one_twelfth * (cells[0] + neg_five * cells[1] + thirteen * cells[2] + three * cells[3])
    }

    /// 4th-order ENO stencil 1 (right): {u_{j-2}, u_{j-1}, u_j, u_{j+1}}
    #[must_use]
    fn eno4_stencil_1_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let seven = from_f64::<T>(7.0);

        one_twelfth * (-cells[1] + seven * cells[2] + seven * cells[3] - cells[4])
    }

    /// 4th-order ENO stencil 2 (right): {u_{j-1}, u_j, u_{j+1}, u_{j+2}}
    #[must_use]
    fn eno4_stencil_2_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let neg_five = from_f64::<T>(-5.0);
        let thirteen = from_f64::<T>(13.0);
        let three = from_f64::<T>(3.0);

        one_twelfth * (three * cells[2] + thirteen * cells[3] + neg_five * cells[4] + cells[5])
    }

    /// 4th-order ENO stencil 3 (right): {u_j, u_{j+1}, u_{j+2}, u_{j+3}}
    #[must_use]
    fn eno4_stencil_3_right(&self, cells: &[T; 7]) -> T {
        let one_twelfth = from_f64::<T>(1.0 / 12.0);
        let neg_three = from_f64::<T>(-3.0);
        let thirteen = from_f64::<T>(13.0);
        let neg_twenty_three = from_f64::<T>(-23.0);
        let twenty_five = from_f64::<T>(25.0);

        one_twelfth
            * (twenty_five * cells[3]
                + neg_twenty_three * cells[4]
                + thirteen * cells[5]
                + neg_three * cells[6])
    }

    /// Smoothness indicator beta0 for WENO7
    #[must_use]
    fn smoothness_indicator_0(&self, cells: &[T; 7]) -> T {
        let one_240 = from_f64::<T>(1.0 / 240.0);
        let c2107 = from_f64::<T>(2107.0);
        let c9402 = from_f64::<T>(9402.0);
        let c7042 = from_f64::<T>(7042.0);
        let c1854 = from_f64::<T>(1854.0);
        let c11003 = from_f64::<T>(11003.0);
        let c17246 = from_f64::<T>(17246.0);
        let c4642 = from_f64::<T>(4642.0);
        let c7043 = from_f64::<T>(7043.0);
        let c3882 = from_f64::<T>(3882.0);
        let c547 = from_f64::<T>(547.0);

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
        let one_240 = from_f64::<T>(1.0 / 240.0);
        let c547 = from_f64::<T>(547.0);
        let c2522 = from_f64::<T>(2522.0);
        let c1922 = from_f64::<T>(1922.0);
        let c494 = from_f64::<T>(494.0);
        let c3443 = from_f64::<T>(3443.0);
        let c5966 = from_f64::<T>(5966.0);
        let c1602 = from_f64::<T>(1602.0);
        let c2843 = from_f64::<T>(2843.0);
        let c1642 = from_f64::<T>(1642.0);
        let c267 = from_f64::<T>(267.0);

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
        let one_240 = from_f64::<T>(1.0 / 240.0);
        let c267 = from_f64::<T>(267.0);
        let c1642 = from_f64::<T>(1642.0);
        let c1602 = from_f64::<T>(1602.0);
        let c494 = from_f64::<T>(494.0);
        let c2843 = from_f64::<T>(2843.0);
        let c5966 = from_f64::<T>(5966.0);
        let c1922 = from_f64::<T>(1922.0);
        let c3443 = from_f64::<T>(3443.0);
        let c2522 = from_f64::<T>(2522.0);
        let c547 = from_f64::<T>(547.0);

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
        let one_240 = from_f64::<T>(1.0 / 240.0);
        let c547 = from_f64::<T>(547.0);
        let c3882 = from_f64::<T>(3882.0);
        let c4642 = from_f64::<T>(4642.0);
        let c1854 = from_f64::<T>(1854.0);
        let c7043 = from_f64::<T>(7043.0);
        let c17246 = from_f64::<T>(17246.0);
        let c7042 = from_f64::<T>(7042.0);
        let c11003 = from_f64::<T>(11003.0);
        let c9402 = from_f64::<T>(9402.0);
        let c2107 = from_f64::<T>(2107.0);

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

impl<T: RealField + Copy + FloatElement> Default for WENO7<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

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
