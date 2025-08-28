//! Interpolation methods for CFD simulations.
//!
//! This module provides various interpolation algorithms optimized for CFD applications
//! with support for both regular and irregular grids.

mod cubic_spline;
mod lagrange;
mod linear;
mod traits;

pub use cubic_spline::CubicSplineInterpolation;
pub use lagrange::LagrangeInterpolation;
pub use linear::LinearInterpolation;
pub use traits::Interpolation;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::error::Result;

    #[test]
    fn test_linear_interpolation() -> Result<()> {
        let x_data = vec![0.0, 1.0, 2.0];
        let y_data = vec![0.0, 1.0, 4.0];

        let interp = LinearInterpolation::new(x_data, y_data)?;

        // Test exact points
        assert_eq!(interp.interpolate(0.0)?, 0.0);
        assert_eq!(interp.interpolate(1.0)?, 1.0);
        assert_eq!(interp.interpolate(2.0)?, 4.0);

        // Test interpolation
        assert_relative_eq!(interp.interpolate(0.5)?, 0.5, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.5)?, 2.5, epsilon = 1e-10);

        Ok(())
    }

    #[test]
    fn test_cubic_spline_interpolation() -> Result<()> {
        let x_data = vec![0.0, 1.0, 2.0, 3.0];
        let y_data = vec![0.0, 1.0, 4.0, 9.0];

        let spline = CubicSplineInterpolation::new(x_data, y_data)?;

        // Test exact points
        assert_relative_eq!(spline.interpolate(0.0)?, 0.0, epsilon = 1e-10);
        assert_relative_eq!(spline.interpolate(2.0)?, 4.0, epsilon = 1e-10);

        // Test smoothness (approximate for quadratic)
        assert_relative_eq!(spline.interpolate(1.5)?, 2.25, epsilon = 0.1);
        assert_relative_eq!(spline.interpolate(2.5)?, 6.25, epsilon = 0.1);

        Ok(())
    }

    #[test]
    fn test_lagrange_interpolation() -> Result<()> {
        let x_data = vec![0.0, 1.0, 2.0];
        let y_data = vec![1.0, 3.0, 7.0];

        let interp = LagrangeInterpolation::new(x_data, y_data)?;

        // Test exact points
        assert_relative_eq!(interp.interpolate(0.0)?, 1.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(1.0)?, 3.0, epsilon = 1e-10);
        assert_relative_eq!(interp.interpolate(2.0)?, 7.0, epsilon = 1e-10);

        // Test interpolation (quadratic: y = x^2 + x + 1)
        assert_relative_eq!(interp.interpolate(0.5)?, 1.75, epsilon = 1e-10);

        Ok(())
    }
}
