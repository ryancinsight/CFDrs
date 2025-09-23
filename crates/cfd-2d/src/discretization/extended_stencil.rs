//! Extended stencil schemes for high-order discretization
//!
//! Provides proper implementation of schemes requiring more than nearest neighbors

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Extended stencil trait for schemes needing wider access
pub trait ExtendedStencilScheme<T: RealField + Copy> {
    /// Compute face value using extended stencil
    ///
    /// # Arguments
    /// * `values` - Array of cell values [`φ_UU`, `φ_U`, `φ_C`, `φ_D`, `φ_DD`]
    /// * `velocity` - Face velocity (for upwind determination)
    /// * `positions` - Cell positions for non-uniform grids
    fn face_value(&self, values: &[T; 5], velocity: T, positions: Option<&[T; 5]>) -> T;
}

/// QUICK scheme with proper 3-point upstream stencil
///
/// Reference: Leonard, B.P. (1979). Computer Methods in Applied Mechanics and Engineering, 19(1), 59-98
pub struct QuickScheme;

impl<T: RealField + Copy + FromPrimitive> ExtendedStencilScheme<T> for QuickScheme {
    fn face_value(&self, values: &[T; 5], velocity: T, _positions: Option<&[T; 5]>) -> T {
        // Constants from Leonard (1979)
        const SIX_EIGHTHS: f64 = 6.0 / 8.0;
        const THREE_EIGHTHS: f64 = 3.0 / 8.0;
        const ONE_EIGHTH: f64 = 1.0 / 8.0;

        let six_eighths = T::from_f64(SIX_EIGHTHS).unwrap_or_else(T::one);
        let three_eighths = T::from_f64(THREE_EIGHTHS).unwrap_or_else(T::zero);
        let one_eighth = T::from_f64(ONE_EIGHTH).unwrap_or_else(T::zero);

        if velocity > T::zero() {
            // Flow from left to right: use φ_U, φ_C, φ_D
            // φ_f = 6/8 * φ_C + 3/8 * φ_D - 1/8 * φ_U
            six_eighths * values[2] + three_eighths * values[3] - one_eighth * values[1]
        } else {
            // Flow from right to left: use φ_D, φ_C, φ_U
            // φ_f = 6/8 * φ_C + 3/8 * φ_U - 1/8 * φ_D
            six_eighths * values[2] + three_eighths * values[1] - one_eighth * values[3]
        }
    }
}

/// MUSCL scheme with van Leer limiter
///
/// Reference: van Leer, B. (1979). Journal of Computational Physics, 32(1), 101-136
pub struct MusclScheme;

impl<T: RealField + Copy + FromPrimitive> ExtendedStencilScheme<T> for MusclScheme {
    fn face_value(&self, values: &[T; 5], velocity: T, _positions: Option<&[T; 5]>) -> T {
        let half = T::from_f64(0.5).unwrap_or_else(T::zero);

        // Compute gradients
        let grad_left = values[2] - values[1];
        let grad_right = values[3] - values[2];

        // Van Leer limiter
        let limiter = |r: T| -> T {
            if r <= T::zero() {
                T::zero()
            } else {
                let two = T::from_f64(2.0).unwrap_or_else(T::one);
                (two * r) / (T::one() + r)
            }
        };

        if velocity > T::zero() {
            let r = grad_left / (grad_right + T::default_epsilon());
            values[2] + half * limiter(r) * grad_right
        } else {
            let r = grad_right / (grad_left + T::default_epsilon());
            values[2] - half * limiter(r) * grad_left
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_scheme() {
        let scheme = QuickScheme;
        let values = [1.0, 2.0, 3.0, 4.0, 5.0];

        // Test positive velocity
        let result = scheme.face_value(&values, 1.0, None);
        let expected: f64 = 6.0 / 8.0 * 3.0 + 3.0 / 8.0 * 4.0 - 1.0 / 8.0 * 2.0;
        assert!((result - expected).abs() < 1e-10f64);

        // Test negative velocity
        let result = scheme.face_value(&values, -1.0, None);
        let expected: f64 = 6.0 / 8.0 * 3.0 + 3.0 / 8.0 * 2.0 - 1.0 / 8.0 * 4.0;
        assert!((result - expected).abs() < 1e-10f64);
    }
}
