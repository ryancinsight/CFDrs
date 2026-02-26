//! Interpolation methods for CFD simulations.
//!
//! This module provides various interpolation algorithms optimized for CFD applications
//! with support for both regular and irregular grids.
//!
//! ## Theorem — Runge Phenomenon (Equidistant Node Instability)
//!
//! **Theorem (Runge 1901)**: For the function f(x) = 1/(1 + 25x²) on [-1, 1],
//! polynomial interpolation at N+1 equidistant nodes diverges as N → ∞:
//!
//! ```text
//! max_{x∈[-1,1]} |f(x) - pₙ(x)| → ∞
//! ```
//!
//! **Implication for CFD**: Equidistant high-order polynomial interpolation is
//! numerically unsafe. All high-order interpolation in this module either uses
//! Chebyshev nodes (which eliminate Runge) or is bounded-degree (linear, cubic).
//!
//! ## Theorem — Lebesgue Constant Bound (Chebyshev Nodes)
//!
//! **Theorem**: For Chebyshev nodes of the first kind xᵢ = cos((2i+1)π/(2N+2)),
//! the Lebesgue constant satisfies:
//!
//! ```text
//! Λₙ = ‖Σᵢ |ℓᵢ(x)|‖_∞ ≤ (2/π) ln(N+1) + 1
//! ```
//!
//! where ℓᵢ(x) are the Lagrange basis polynomials. This O(log N) growth (vs
//! the exponential growth for equidistant nodes) guarantees near-optimal
//! approximation in the L∞ sense.
//!
//! ## Theorem — Cubic Spline Smoothness
//!
//! **Theorem**: The natural cubic spline interpolant S(x) minimises the bending energy
//!
//! ```text
//! ∫_a^b [S''(x)]² dx   subject to S(xᵢ) = yᵢ  ∀i
//! ```
//!
//! among all functions in C²([a,b]). This is the variational characterisation of
//! the cubic spline as the "minimum curvature" interpolant.
//!
//! ## Theorem — Error Bound for Cubic Spline
//!
//! **Theorem**: For f ∈ C⁴([a,b]) with mesh spacing h = maxᵢ(xᵢ₊₁ - xᵢ):
//!
//! ```text
//! ‖f - S‖_∞ ≤ (5/384) h⁴ ‖f⁽⁴⁾‖_∞
//! ```
//!
//! ## Invariants
//! - Interpolation nodes must be strictly increasing: x₀ < x₁ < … < xₙ.
//! - Evaluation outside [x₀, xₙ] is extrapolation — the module returns an error.
//! - For linear interpolation, the result is exact at the data points.

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

    // --- Adversarial / Robustness Tests ------------------------------------------

    /// Verify extrapolation returns Err, not a silently wrong value.
    #[test]
    fn test_linear_extrapolation_rejected() {
        let x_data = vec![0.0, 1.0, 2.0];
        let y_data = vec![0.0, 1.0, 4.0];
        let interp = LinearInterpolation::new(x_data, y_data)
            .expect("construction must succeed for valid data");
        // Query below domain minimum
        let result = interp.interpolate(-0.1);
        assert!(result.is_err(), "extrapolation below domain must return Err");
        // Query above domain maximum
        let result = interp.interpolate(2.1);
        assert!(result.is_err(), "extrapolation above domain must return Err");
    }

    /// Verify that a single-point dataset (degenerate case) is handled.
    #[test]
    fn test_linear_single_point_rejected() {
        let x_data = vec![1.0];
        let y_data = vec![3.0];
        let result = LinearInterpolation::new(x_data, y_data);
        // A single point cannot define an interpolation scheme — must fail gracefully.
        assert!(result.is_err(), "single-point interpolation must return Err");
    }

    /// Verify non-strictly-increasing nodes are rejected.
    #[test]
    fn test_interpolation_duplicate_nodes_rejected() {
        let x_data = vec![0.0, 1.0, 1.0, 2.0]; // duplicate node at 1.0
        let y_data = vec![0.0, 1.0, 2.0, 3.0];
        let result = LinearInterpolation::new(x_data, y_data);
        assert!(result.is_err(), "duplicate nodes must return Err");
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

    /// Natural cubic spline on quadratic data: not exact because natural BCs
    /// (S''(x₀) = S''(xₙ) = 0) conflict with f''(x) = 2 for f(x) = x².
    /// Natural splines reproduce only polynomials of degree ≤ 1 exactly.
    /// Still converges: ‖f − S‖_∞ ≤ (5/384)h⁴‖f⁽⁴⁾‖_∞ = 0 for f(x)=x²
    /// at interior knots, with O(h²) boundary pollution.
    #[test]
    fn test_cubic_spline_quadratic_exactness() -> Result<()> {
        // f(x) = x² — natural spline approximates well but is not exact
        // due to S''(0) = 0 ≠ f''(0) = 2 boundary mismatch.
        let xs = vec![0.0_f64, 0.5, 1.0, 1.5, 2.0];
        let ys: Vec<f64> = xs.iter().map(|&x| x * x).collect();
        let spline = CubicSplineInterpolation::new(xs.clone(), ys)?;
        // Interior points are better approximated than near-boundary points.
        for xi in [0.25, 0.75, 1.25, 1.75_f64] {
            let exact = xi * xi;
            let approx = spline.interpolate(xi)?;
            assert_relative_eq!(approx, exact, epsilon = 0.05);
        }
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
