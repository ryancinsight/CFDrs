//! Finite difference operators for numerical differentiation.

use super::schemes::FiniteDifferenceScheme;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Finite difference operator for 1D problems
pub struct FiniteDifference<T: RealField + Copy> {
    scheme: FiniteDifferenceScheme,
    spacing: T,
}

impl<T: RealField + FloatElement + Copy> FiniteDifference<T> {
    /// Create a finite difference operator
    pub fn new(scheme: FiniteDifferenceScheme, spacing: T) -> Self {
        Self { scheme, spacing }
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
    ///
    /// # Errors
    /// Returns an error if the input array has fewer than 2 points
    pub fn first_derivative(&self, values: &[T]) -> Result<Array1<T>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string(),
            ));
        }

        let n = values.len();
        let mut result = vector_zeros(n);
        let inv_spacing = <T as NumericElement>::ONE / self.spacing;

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
                let two_inv_spacing = inv_spacing / from_f64::<T>(2.0);
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

                let two = from_f64::<T>(2.0);
                let three = from_f64::<T>(3.0);
                let four = from_f64::<T>(4.0);

                // Use forward difference for first n-2 points
                for i in 0..n.saturating_sub(2) {
                    result[i] = (-three * values[i] + four * values[i + 1] - values[i + 2])
                        / (two * self.spacing);
                }

                // Use central difference for remaining points
                for i in n.saturating_sub(2)..n {
                    if i > 0 && i < n - 1 {
                        result[i] = (values[i + 1] - values[i - 1]) / (two * self.spacing);
                    }
                }
            }
            FiniteDifferenceScheme::BackwardSecondOrder => {
                if n < 3 {
                    return Err(Error::InvalidConfiguration(
                        "Need at least 3 points for second-order backward difference".to_string(),
                    ));
                }

                let two = from_f64::<T>(2.0);
                let three = from_f64::<T>(3.0);
                let four = from_f64::<T>(4.0);

                // Use central difference for first points
                for i in 0..2 {
                    if i > 0 && i < n - 1 {
                        result[i] = (values[i + 1] - values[i - 1]) / (two * self.spacing);
                    }
                }

                // Use backward difference for remaining points
                for i in 2..n {
                    result[i] = (values[i - 2] - four * values[i - 1] + three * values[i])
                        / (two * self.spacing);
                }
            }
        }

        Ok(result)
    }

    /// Compute second derivative using central differences
    ///
    /// # Errors
    /// Returns an error if the input array has fewer than 3 points
    pub fn second_derivative(&self, values: &[T]) -> Result<Array1<T>> {
        if values.len() < 3 {
            return Err(Error::InvalidConfiguration(
                "Need at least 3 points for second derivative".to_string(),
            ));
        }

        let n = values.len();
        let mut result = vector_zeros(n);
        let h_squared = self.spacing * self.spacing;
        let two = from_f64::<T>(2.0);

        // Use forward difference for first point
        result[0] = (values[2] - two * values[1] + values[0]) / h_squared;

        // Central difference for interior points
        for i in 1..n - 1 {
            result[i] = (values[i + 1] - two * values[i] + values[i - 1]) / h_squared;
        }

        // Use backward difference for last point
        result[n - 1] = (values[n - 1] - two * values[n - 2] + values[n - 3]) / h_squared;

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

impl<T: RealField + FloatElement + Copy> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::central(<T as NumericElement>::ONE)
    }
}

impl FiniteDifference<f32> {
    /// Compute first derivative using the f32 SIMD-friendly path.
    ///
    /// # Errors
    /// Returns an error if the input array has fewer than 2 points.
    pub fn first_derivative_simd(&self, values: &[f32]) -> Result<Vec<f32>> {
        if values.len() < 2 {
            return Err(Error::InvalidConfiguration(
                "Need at least 2 points for differentiation".to_string(),
            ));
        }

        let n = values.len();
        let mut result = vec![0.0f32; n];
        let inv_spacing = 1.0f32 / self.spacing;

        if self.scheme == FiniteDifferenceScheme::Central {
            if n > 2 {
                let scale = inv_spacing * 0.5;
                for i in 1..n - 1 {
                    result[i] = (values[i + 1] - values[i - 1]) * scale;
                }
            }

            result[0] = (values[1] - values[0]) * inv_spacing;
            result[n - 1] = (values[n - 1] - values[n - 2]) * inv_spacing;
        } else {
            let scalar_result = self.first_derivative(values)?;
            for (i, val) in scalar_result.iter().enumerate() {
                result[i] = *val;
            }
        }

        Ok(result)
    }
}

/// Compute 1D differentiation using central differences
///
/// # Errors
/// Returns an error if the input array has fewer than 2 points
pub fn differentiate_1d<T: RealField + FloatElement + Copy>(
    values: &[T],
    spacing: T,
) -> Result<Array1<T>> {
    FiniteDifference::central(spacing).first_derivative(values)
}

/// Compute 2D differentiation (gradient) using central differences
///
/// # Errors
/// Returns an error if field dimensions are invalid or insufficient data points
pub fn differentiate_2d<T: RealField + FloatElement + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Result<(Vec<T>, Vec<T>)> {
    use crate::differentiation::Gradient;

    let grad = Gradient::new(dx, dy, <T as NumericElement>::ONE);
    let gradients = grad.gradient_2d(field, nx, ny)?;

    let grad_x: Vec<T> = gradients.iter().map(|g| g.x).collect();
    let grad_y: Vec<T> = gradients.iter().map(|g| g.y).collect();

    Ok((grad_x, grad_y))
}

/// Compute 2D Laplacian using central differences
///
/// # Errors
/// Returns an error if field dimensions are invalid or insufficient data points
#[allow(clippy::similar_names)] // Mathematical derivatives use standard notation
pub fn laplacian_2d<T: RealField + FloatElement + Copy>(
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

    let mut laplacian = vec![<T as NumericElement>::ZERO; nx * ny];
    let dx2 = dx * dx;
    let dy2 = dy * dy;
    let two = from_f64::<T>(2.0);

    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;

            // ∂²f/∂x²
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;

            // ∂²f/∂y²
            let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;

            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Handle boundaries with forward/backward differences
    // Top and bottom boundaries
    for i in 0..nx {
        // Bottom boundary (j=0)
        let idx = i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;
            let d2fdy2 = (field[idx + 2 * nx] - two * field[idx + nx] + field[idx]) / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }

        // Top boundary (j=ny-1)
        let idx = (ny - 1) * nx + i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;
            let d2fdy2 = (field[idx] - two * field[idx - nx] + field[idx - 2 * nx]) / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Left and right boundaries
    for j in 1..ny - 1 {
        // Left boundary (i=0)
        let idx = j * nx;
        let d2fdx2 = (field[idx + 2] - two * field[idx + 1] + field[idx]) / dx2;
        let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;

        // Right boundary (i=nx-1)
        let idx = j * nx + (nx - 1);
        let d2fdx2 = (field[idx] - two * field[idx - 1] + field[idx - 2]) / dx2;
        let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;
    }

    Ok(laplacian)
}

fn vector_zeros<T: NumericElement>(len: usize) -> Array1<T> {
    Array1::from_elem([len], T::ZERO)
}
