//! Rhie-Chow momentum interpolation
//!
//! Reference: Rhie, C.M. and Chow, W.L. (1983). "Numerical study of the turbulent
//! flow past an airfoil with trailing edge separation". AIAA Journal, 21(11), 1525-1532.

use eunomia::FloatElement;
use eunomia::RealField;

/// Rhie-Chow interpolation for pressure-velocity coupling
///
/// This interpolation prevents checkerboard pressure oscillations in collocated grids
/// by adding a pressure-gradient-based correction to the interpolated velocity.
pub struct RhieChowInterpolation<T: RealField + Copy> {
    /// Grid spacing in x-direction
    dx: T,
    /// Grid spacing in y-direction
    dy: T,
    /// Relaxation factor
    alpha: T,
}

impl<T: RealField + FloatElement + Copy> RhieChowInterpolation<T> {
    /// Create new Rhie-Chow interpolator
    pub fn new(dx: T, dy: T) -> Self {
        Self {
            dx,
            dy,
            alpha: <T as FloatElement>::from_f64(0.8),
        }
    }

    /// Helper to get the constant 2.0 safely
    #[inline]
    fn two() -> T {
        <T as FloatElement>::from_f64(2.0)
    }

    /// Interpolate velocity to face with pressure correction
    ///
    /// `u_f` = `ū_f` - `D_f` * (∇p)_f
    ///
    /// where:
    /// - `u_f` is the face velocity
    /// - `ū_f` is the linearly interpolated velocity
    /// - `D_f` is the face diffusion coefficient
    /// - (∇p)_f is the pressure gradient at the face
    pub fn interpolate_u_face(
        &self,
        u_p: T, // velocity at cell P
        u_e: T, // velocity at cell E
        p_p: T, // pressure at cell P
        p_e: T, // pressure at cell E
        a_p: T, // diagonal coefficient for cell P
        a_e: T, // diagonal coefficient for cell E
        dt: T,  // time step
    ) -> T {
        // Linear interpolation of velocity
        let two = Self::two();
        let u_bar = (u_p + u_e) / two;

        // Face diffusion coefficient (harmonic mean)
        let d_f = dt * two / (a_p + a_e);

        // Pressure gradient at face
        let dp_dx = (p_e - p_p) / self.dx;

        // Rhie-Chow correction with relaxation
        let correction = d_f * dp_dx;
        u_bar - self.alpha * correction
    }

    /// Interpolate v-velocity to face with pressure correction
    pub fn interpolate_v_face(
        &self,
        v_p: T, // velocity at cell P
        v_n: T, // velocity at cell N
        p_p: T, // pressure at cell P
        p_n: T, // pressure at cell N
        a_p: T, // diagonal coefficient for cell P
        a_n: T, // diagonal coefficient for cell N
        dt: T,  // time step
    ) -> T {
        // Linear interpolation of velocity
        let two = Self::two();
        let v_bar = (v_p + v_n) / two;

        // Face diffusion coefficient
        let d_f = dt * two / (a_p + a_n);

        // Pressure gradient at face
        let dp_dy = (p_n - p_p) / self.dy;

        // Rhie-Chow correction with relaxation
        let correction = d_f * dp_dy;
        v_bar - self.alpha * correction
    }
}

#[cfg(test)]
mod tests {
    use super::RhieChowInterpolation;

    fn assert_close(actual: f64, expected: f64) {
        let scale = actual.abs().max(expected.abs()).max(1.0);
        // Straight-line interpolation has fewer than eight rounded operations.
        let tolerance = f64::EPSILON * scale * 8.0;
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual}, expected={expected}, tolerance={tolerance}"
        );
    }

    #[test]
    fn u_face_interpolation_applies_pressure_correction() {
        let interpolation = RhieChowInterpolation::new(0.5_f64, 0.25);

        let actual = interpolation.interpolate_u_face(
            2.0,  // u_p
            4.0,  // u_e
            10.0, // p_p
            14.0, // p_e
            3.0,  // a_p
            5.0,  // a_e
            0.2,  // dt
        );

        let u_bar = f64::midpoint(2.0, 4.0);
        let d_f = 0.2 * 2.0 / (3.0 + 5.0);
        let dp_dx = (14.0 - 10.0) / 0.5;
        let expected = u_bar - 0.8 * d_f * dp_dx;
        assert_close(actual, expected);
    }

    #[test]
    fn v_face_interpolation_applies_pressure_correction() {
        let interpolation = RhieChowInterpolation::new(0.5_f64, 0.25);

        let actual = interpolation.interpolate_v_face(
            -1.0, // v_p
            3.0,  // v_n
            8.0,  // p_p
            9.0,  // p_n
            2.0,  // a_p
            6.0,  // a_n
            0.4,  // dt
        );

        let v_bar = f64::midpoint(-1.0, 3.0);
        let d_f = 0.4 * 2.0 / (2.0 + 6.0);
        let dp_dy = (9.0 - 8.0) / 0.25;
        let expected = v_bar - 0.8 * d_f * dp_dy;
        assert_close(actual, expected);
    }
}
