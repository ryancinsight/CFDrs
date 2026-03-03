//! Error estimation methods for adaptive mesh refinement.
//!
//! Provides Richardson extrapolation, adjoint-based, residual-based,
//! and smoothness-based error estimators.

use nalgebra::DMatrix;

use super::AdaptiveMeshRefinement;

impl AdaptiveMeshRefinement {
    /// Compute sophisticated error estimate using multiple methods
    pub(super) fn compute_sophisticated_error_estimate(
        solution: &DMatrix<f64>,
        i: usize,
        j: usize,
        threshold: f64,
    ) -> f64 {
        // Combine multiple error estimation techniques
        let richardson_error = Self::compute_richardson_error_estimate(solution, i, j);
        let adjoint_error = Self::compute_adjoint_error_estimate(solution, i, j);
        let residual_error = Self::compute_residual_error_estimate(solution, i, j);
        let smoothness_error = Self::compute_smoothness_error_estimate(solution, i, j);

        // Weighted combination of error estimates
        let combined_error = 0.3 * richardson_error
            + 0.3 * adjoint_error
            + 0.2 * residual_error
            + 0.2 * smoothness_error;

        // Apply threshold normalization
        if combined_error > threshold {
            combined_error / threshold
        } else {
            0.0
        }
    }

    /// Richardson extrapolation-based error estimate
    fn compute_richardson_error_estimate(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        // Need at least 2x2 stencil for Richardson extrapolation
        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Fine grid approximation (current solution)
        let fine = solution[(i, j)];

        // Coarse grid approximation (average of neighbors)
        let coarse = 0.25
            * (solution[(i - 1, j)]
                + solution[(i + 1, j)]
                + solution[(i, j - 1)]
                + solution[(i, j + 1)]);

        // Richardson error estimate (assuming second-order accuracy)
        let error = (fine - coarse).abs() / 3.0;

        // Scale by local solution magnitude for relative error
        if fine.abs() > 1e-10 {
            error / fine.abs()
        } else {
            error
        }
    }

    /// Adjoint-based error estimate (simplified)
    fn compute_adjoint_error_estimate(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute local Hessian (second derivatives)
        let d2u_dx2 = solution[(i + 1, j)] - 2.0 * solution[(i, j)] + solution[(i - 1, j)];
        let d2u_dy2 = solution[(i, j + 1)] - 2.0 * solution[(i, j)] + solution[(i, j - 1)];
        let d2u_dxdy = 0.25
            * (solution[(i + 1, j + 1)] - solution[(i + 1, j - 1)] - solution[(i - 1, j + 1)]
                + solution[(i - 1, j - 1)]);

        // Adjoint weight (simplified - based on solution curvature)
        let adjoint_weight =
            (d2u_dx2 * d2u_dx2 + d2u_dy2 * d2u_dy2 + 2.0 * d2u_dxdy * d2u_dxdy).sqrt();

        // Error estimate based on local truncation error weighted by adjoint
        let local_truncation_error = (d2u_dx2.abs() + d2u_dy2.abs()) / 12.0;

        adjoint_weight * local_truncation_error
    }

    /// Residual-based error estimate
    fn compute_residual_error_estimate(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute discrete Laplacian (residual of Poisson equation)
        let laplacian = solution[(i + 1, j)]
            + solution[(i - 1, j)]
            + solution[(i, j + 1)]
            + solution[(i, j - 1)]
            - 4.0 * solution[(i, j)];

        // Residual magnitude
        laplacian.abs()
    }

    /// Smoothness-based error estimate
    fn compute_smoothness_error_estimate(solution: &DMatrix<f64>, i: usize, j: usize) -> f64 {
        let (nx, ny) = solution.shape();

        if i == 0 || i >= nx - 1 || j == 0 || j >= ny - 1 {
            return 0.0;
        }

        // Compute solution variation in different directions
        let dx_plus = solution[(i + 1, j)] - solution[(i, j)];
        let dx_minus = solution[(i, j)] - solution[(i - 1, j)];
        let dy_plus = solution[(i, j + 1)] - solution[(i, j)];
        let dy_minus = solution[(i, j)] - solution[(i, j - 1)];

        // Smoothness indicator (based on total variation diminishing)
        let smoothness_x = if dx_plus.abs() + dx_minus.abs() > 1e-10 {
            (dx_plus - dx_minus).abs() / (dx_plus.abs() + dx_minus.abs())
        } else {
            0.0
        };

        let smoothness_y = if dy_plus.abs() + dy_minus.abs() > 1e-10 {
            (dy_plus - dy_minus).abs() / (dy_plus.abs() + dy_minus.abs())
        } else {
            0.0
        };

        // Combined smoothness error
        f64::midpoint(smoothness_x, smoothness_y)
    }
}
