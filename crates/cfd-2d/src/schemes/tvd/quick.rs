//! QUICK (Quadratic Upstream Interpolation for Convective Kinematics) scheme.

use nalgebra::RealField;

/// QUICK scheme (Leonard, 1979) — third-order accurate interpolation
/// based on quadratic upstream interpolation. Combines second-order central
/// differencing with third-order upwind-biased interpolation.
pub struct QUICKScheme<T: RealField + Copy> {
    /// Courant number for stability analysis
    courant_max: T,
}

impl<T: RealField + Copy> Default for QUICKScheme<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy> QUICKScheme<T> {
    /// Create new QUICK scheme with default Courant limit of 0.8
    #[must_use]
    pub fn new() -> Self {
        Self {
            courant_max: crate::schemes::constants::to_realfield::<T>(
                crate::schemes::constants::CFL_QUICK_SCHEME,
            ),
        }
    }

    /// Compute face value using QUICK interpolation
    ///
    /// For uniform grid with flow from left to right:
    /// `φ_face` = 6/8 * `φ_U` + 3/8 * `φ_C` - 1/8 * `φ_UU`
    /// where U = upstream, C = central, UU = far upstream
    pub fn interpolate_face(
        &self,
        phi_uu: T, // Far upstream value
        phi_u: T,  // Upstream value
        phi_c: T,  // Central value
        phi_d: T,  // Downstream value
        velocity: T,
    ) -> T {
        if velocity > T::zero() {
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_u + three_eighths * phi_c - one_eighth * phi_uu
        } else {
            let six_eighths = T::from_f64(0.75).unwrap_or_else(T::one);
            let three_eighths = T::from_f64(0.375).unwrap_or_else(T::one);
            let one_eighth = T::from_f64(0.125).unwrap_or_else(T::one);

            six_eighths * phi_c + three_eighths * phi_u - one_eighth * phi_d
        }
    }

    /// Check Courant number for stability
    pub fn is_stable(&self, courant: T) -> bool {
        courant <= self.courant_max
    }
}
