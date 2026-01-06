//! Helper functions for Spalart-Allmaras model
//!
//! Extracted utility functions to maintain <500 line limit per file

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

    if x.abs() < T::from_f64(EPSILON_MIN).unwrap_or_else(T::one) {
        return T::zero();
    }

    // Initial guess: x^(1/3) ≈ x/3
    let three = T::from_f64(THREE).unwrap_or_else(T::one);
    let two = T::from_f64(TWO).unwrap_or_else(T::one);

    let mut guess = x / three;

    // Newton-Raphson: y_{n+1} = (2*y_n + x/y_n²) / 3
    for _ in 0..10 {
        let guess_sq = guess * guess;
        let new_guess = (two * guess + x / guess_sq) / three;
        if (new_guess - guess).abs() < T::from_f64(1e-10).unwrap_or_else(T::one) {
            break;
        }
        guess = new_guess;
    }

    guess
}

/// Calculate wall distance field for 2D rectangular domain
///
/// Simplified 2D wall distance: minimum distance to domain boundaries
/// In production code, use Eikonal equation solver for complex geometries
#[instrument(skip(dx, dy))]
pub fn wall_distance_field_2d<T: RealField + Copy + FromPrimitive>(
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Vec<T> {
    const EPSILON_MIN: f64 = 1e-10;
    let mut distances = Vec::with_capacity(nx * ny);

    for j in 0..ny {
        for i in 0..nx {
            // Distance from boundaries (assuming rectangular domain)
            let i_t = T::from_usize(i).unwrap_or_else(T::zero);
            let j_t = T::from_usize(j).unwrap_or_else(T::zero);
            let nx_t = T::from_usize(nx - 1).unwrap_or_else(T::one);
            let ny_t = T::from_usize(ny - 1).unwrap_or_else(T::one);

            let dist_west = i_t * dx;
            let dist_east = (nx_t - i_t) * dx;
            let dist_south = j_t * dy;
            let dist_north = (ny_t - j_t) * dy;

            // Minimum distance to any wall
            let min_dist = dist_west
                .min(dist_east)
                .min(dist_south)
                .min(dist_north)
                .max(T::from_f64(EPSILON_MIN).unwrap_or_else(T::one));

            distances.push(min_dist);
        }
    }

    distances
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
        assert_relative_eq!(distances[corner_idx], 0.0, epsilon = 1e-10);
    }
}
