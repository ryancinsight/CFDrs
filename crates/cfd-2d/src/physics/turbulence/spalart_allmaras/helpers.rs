//! Helper functions for Spalart-Allmaras model.
//!
//! Extracted utility functions to maintain <500 line limit per file and to
//! share the rectangular wall-distance field with the turbulence boundary manager.
//!
//! # Theorem
//! On an orthogonal rectangular grid with cell-centered storage, the shared wall-distance field is
//! non-negative and equals the minimum center-to-wall distance for each cell.
//!
//! **Proof sketch**:
//! [`wall_distance_field_2d`] delegates to the turbulence boundary manager, which evaluates the
//! minimum of the four center-to-boundary distances
//! $x_c$, $L_x - x_c$, $y_c$, and $L_y - y_c$ for each cell center $(x_c, y_c)$. Each term is
//! non-negative on the closed domain, so the minimum is also non-negative. Boundary-adjacent cells
//! therefore evaluate to half a cell spacing in the wall-normal direction, and interior cells grow
//! monotonically with distance from the nearest boundary on a uniform rectangle.

use crate::physics::turbulence::boundary_conditions::TurbulenceBoundaryManager;
use cfd_core::physics::constants::mathematical::numeric::{THREE, TWO};
use nalgebra::RealField;
use num_traits::FromPrimitive;
use tracing::instrument;

/// Compute cube root using Newton-Raphson iteration
///
/// Needed because `RealField` trait doesn't provide cbrt directly
#[instrument(skip(x))]
pub fn cbrt<T: RealField + Copy + FromPrimitive>(x: T) -> T {
    const EPSILON_MIN: f64 = 1e-10;

    if x.abs() < T::from_f64(EPSILON_MIN).expect("analytical constant conversion") {
        return T::zero();
    }

    // Initial guess: x^(1/3) ≈ x/3
    let three = T::from_f64(THREE).expect("analytical constant conversion");
    let two = T::from_f64(TWO).expect("analytical constant conversion");

    let mut guess = x / three;

    // Newton-Raphson: y_{n+1} = (2*y_n + x/y_n²) / 3
    for _ in 0..10 {
        let guess_sq = guess * guess;
        let new_guess = (two * guess + x / guess_sq) / three;
        if (new_guess - guess).abs() < T::from_f64(1e-10).expect("analytical constant conversion") {
            break;
        }
        guess = new_guess;
    }

    guess
}

/// Calculate the rectangular wall-distance field shared with the turbulence
/// boundary manager.
#[instrument(skip(dx, dy))]
pub fn wall_distance_field_2d<T: RealField + Copy + FromPrimitive>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Vec<T> {
    TurbulenceBoundaryManager::new(nx, ny, dx, dy).wall_distances
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_cbrt_function() {
        // Test cube root computation
        assert_relative_eq!(cbrt(8.0_f64), 2.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(27.0_f64), 3.0, epsilon = 1e-8);
        assert_relative_eq!(cbrt(1.0_f64), 1.0, epsilon = 1e-8);
        assert_eq!(cbrt(0.0_f64), 0.0);
    }

    #[test]
    fn test_wall_distance_field_2d() {
        let distances = wall_distance_field_2d(5, 5, 0.1, 0.1);

        assert_eq!(distances.len(), 25);

        // Corner should have minimum distance
        // Center should have maximum distance
        let center_idx = 2 * 5 + 2;
        let corner_idx = 0;

        assert!(distances[center_idx] > distances[corner_idx]);
        assert_relative_eq!(distances[corner_idx], 0.05, epsilon = 1e-10);
    }

    #[test]
    fn test_wall_distance_field_matches_boundary_manager() {
        let distances = wall_distance_field_2d(4, 3, 0.2_f64, 0.1_f64);
        let manager = TurbulenceBoundaryManager::new(4, 3, 0.2_f64, 0.1_f64);

        assert_eq!(distances, manager.wall_distances);
    }
}
