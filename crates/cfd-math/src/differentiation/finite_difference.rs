//! Finite difference operators for numerical differentiation.

use super::schemes::FiniteDifferenceScheme;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// Finite difference operator for 1D problems
pub struct FiniteDifference<T: RealField + Copy> {
    scheme: FiniteDifferenceScheme,
    spacing: T,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> FiniteDifference<T> {
    /// Create a finite difference operator
    pub fn new(scheme: FiniteDifferenceScheme, spacing: T) -> Self {
        Self {
            scheme,
            spacing,
        }
    }

    /// Create central difference operator
    pub fn central(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Central, spacing)
    }

    /// Create forward difference operator
    pub fn forward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Forward, spacing)
    }

    /// Create backward difference operator
    pub fn backward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Backward, spacing)
    }

    /// Compute first derivative using iterator combinators and zero-copy operations
    pub fn first_derivative(&self, values: &[T]) -> Result<DVector<T>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string(),
            ));
        }

        let n = values.len();
        let mut result = DVector::zeros(n);
        let inv_spacing = T::one() / self.spacing;

        match self.scheme {
            FiniteDifferenceScheme::Forward => {
                // Use windows() for efficient forward differences
                values.windows(2).enumerate().for_each(|(i, window)| {
                    result[i] = (window[1] - window[0]) * inv_spacing;
                });

                // Backward difference for last point
                if n > 1 {
                    result[n - 1] = (values[n - 1] - values[n - 2]) * inv_spacing;
                }
            }
            FiniteDifferenceScheme::Backward => {
                // Forward difference for first point
                result[0] = (values[1] - values[0]) * inv_spacing;

                // Use windows() for backward differences: (values[i] - values[i-1])
                values.windows(2).enumerate().for_each(|(i, window)| {
                    result[i + 1] = (window[1] - window[0]) * inv_spacing;
                });
            }
            FiniteDifferenceScheme::Central => {
                // Forward difference for first point
                result[0] = (values[1] - values[0]) * inv_spacing;

                // Central difference using windows(3) for interior points
                let two_inv_spacing = inv_spacing / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                values.windows(3).enumerate().for_each(|(i, window)| {
                    result[i + 1] = (window[2] - window[0]) * two_inv_spacing;
                });

                // Backward difference for last point
                if n > 1 {
                    result[n - 1] = (values[n - 1] - values[n - 2]) * inv_spacing;
                }
            }
            FiniteDifferenceScheme::ForwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order forward difference".to_string(),
                    ));
                }

                let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
                let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
                let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());

                // Use forward difference for first n-2 points
                result
                    .iter_mut()
                    .take(n.saturating_sub(2))
                    .enumerate()
                    .for_each(|(i, r)| {
                        *r = (-three * values[i] + four * values[i + 1] - values[i + 2])
                            / (two * self.spacing);
                    });

                // Use central difference for remaining points
                result
                    .iter_mut()
                    .skip(n.saturating_sub(2))
                    .enumerate()
                    .for_each(|(idx, r)| {
                        let i = idx + n.saturating_sub(2);
                        if i > 0 && i < n - 1 {
                            *r = (values[i + 1] - values[i - 1]) / (two * self.spacing);
                        }
                    });
            }
            FiniteDifferenceScheme::BackwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order backward difference".to_string(),
                    ));
                }

                let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
                let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
                let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());

                // Use central difference for first points
                result.iter_mut().take(2).enumerate().for_each(|(i, r)| {
                    if i > 0 && i < n - 1 {
                        *r = (values[i + 1] - values[i - 1]) / (two * self.spacing);
                    }
                });

                // Use backward difference for remaining points
                result.iter_mut().skip(2).enumerate().for_each(|(idx, r)| {
                    let i = idx + 2;
                    *r = (values[i - 2] - four * values[i - 1] + three * values[i])
                        / (two * self.spacing);
                });
            }
        }

        Ok(result)
    }

    /// Compute first derivative using SIMD acceleration for f32 arrays
    pub fn first_derivative_simd_f32(&self, values: &[f32]) -> Result<Vec<f32>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string(),
            ));
        }

        let n = values.len();
        let mut result = vec![0.0f32; n];
        let inv_spacing = 1.0f32 / (self.spacing.to_subset().unwrap_or(1.0) as f32);

        match self.scheme {
            FiniteDifferenceScheme::Central => {
                // Use SIMD-friendly operations for central differences
                if n > 2 {
                    // Compute differences and scale in one pass
                    let scale = inv_spacing * 0.5;
                    for i in 1..n - 1 {
                        result[i] = (values[i + 1] - values[i - 1]) * scale;
                    }
                }

                // Handle boundaries
                result[0] = (values[1] - values[0]) * inv_spacing;
                result[n - 1] = (values[n - 1] - values[n - 2]) * inv_spacing;
            }
            _ => {
                // Fall back to scalar for other schemes
                let scalar_result = self.first_derivative(
                    &values
                        .iter()
                        .map(|&v| T::from_f32(v).unwrap_or_else(|| T::zero()))
                        .collect::<Vec<_>>(),
                )?;
                for (i, val) in scalar_result.iter().enumerate() {
                    result[i] = val.to_subset().unwrap_or(0.0) as f32;
                }
            }
        }

        Ok(result)
    }

    /// Compute second derivative using central differences
    pub fn second_derivative(&self, values: &[T]) -> Result<DVector<T>> {
        if values.len() < 3 {
            return Err(Error::InvalidConfiguration(
                "Need at least 3 points for second derivative".to_string(),
            ));
        }

        let n = values.len();
        let mut result = DVector::zeros(n);
        let h_squared = self.spacing * self.spacing;

        // Use forward difference for first point
        result[0] = (values[2] - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * values[1]
            + values[0])
            / h_squared;

        // Central difference for interior points
        for i in 1..n - 1 {
            result[i] = (values[i + 1] - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * values[i]
                + values[i - 1])
                / h_squared;
        }

        // Use backward difference for last point
        result[n - 1] = (values[n - 1]
            - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * values[n - 2]
            + values[n - 3])
            / h_squared;

        Ok(result)
    }

    /// Get the finite difference scheme
    pub fn scheme(&self) -> FiniteDifferenceScheme {
        self.scheme
    }

    /// Get the grid spacing
    pub fn spacing(&self) -> T {
        self.spacing
    }
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::central(T::one())
    }
}

/// Compute 1D differentiation using central differences
pub fn differentiate_1d<T: RealField + From<f64> + FromPrimitive + Copy>(
    values: &[T],
    spacing: T,
) -> Result<DVector<T>> {
    FiniteDifference::central(spacing).first_derivative(values)
}

/// Compute 2D differentiation (gradient) using central differences
pub fn differentiate_2d<T: RealField + From<f64> + FromPrimitive + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Result<(Vec<T>, Vec<T>)> {
    use crate::differentiation::Gradient;

    let grad = Gradient::new(dx, dy, T::one());
    let gradients = grad.gradient_2d(field, nx, ny)?;

    let grad_x: Vec<T> = gradients.iter().map(|g| g.x).collect();
    let grad_y: Vec<T> = gradients.iter().map(|g| g.y).collect();

    Ok((grad_x, grad_y))
}

/// Compute 2D Laplacian using central differences
pub fn laplacian_2d<T: RealField + From<f64> + FromPrimitive + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Result<Vec<T>> {
    if field.len() != nx * ny {
        return Err(Error::InvalidConfiguration(
            "Field size doesn't match grid dimensions".to_string(),
        ));
    }

    let mut laplacian = vec![T::zero(); nx * ny];
    let dx2 = dx * dx;
    let dy2 = dy * dy;

    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;

            // ∂²f/∂x²
            let d2fdx2 = (field[idx + 1] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
                + field[idx - 1])
                / dx2;

            // ∂²f/∂y²
            let d2fdy2 = (field[idx + nx] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
                + field[idx - nx])
                / dy2;

            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Handle boundaries with forward/backward differences
    // Top and bottom boundaries
    for i in 0..nx {
        // Bottom boundary (j=0)
        let idx = i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
                + field[idx - 1])
                / dx2;
            let d2fdy2 = (field[idx + 2 * nx]
                - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx + nx]
                + field[idx])
                / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }

        // Top boundary (j=ny-1)
        let idx = (ny - 1) * nx + i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
                + field[idx - 1])
                / dx2;
            let d2fdy2 = (field[idx] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx - nx]
                + field[idx - 2 * nx])
                / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Left and right boundaries
    for j in 1..ny - 1 {
        // Left boundary (i=0)
        let idx = j * nx;
        let d2fdx2 = (field[idx + 2] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx + 1]
            + field[idx])
            / dx2;
        let d2fdy2 = (field[idx + nx] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
            + field[idx - nx])
            / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;

        // Right boundary (i=nx-1)
        let idx = j * nx + (nx - 1);
        let d2fdx2 = (field[idx] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx - 1]
            + field[idx - 2])
            / dx2;
        let d2fdy2 = (field[idx + nx] - T::from_f64(2.0).unwrap_or_else(T::zero) * field[idx]
            + field[idx - nx])
            / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;
    }

    Ok(laplacian)
}
