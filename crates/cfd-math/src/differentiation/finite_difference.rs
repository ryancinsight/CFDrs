//! Finite-difference operators — thin wrapper around the `leto-ops` SSOT.
//!
//! All stencil logic lives in [`leto_ops::FiniteDifference`]; this module
//! re-exports the type alias and adds:
//! - `cfd_core::Error` boundary (via `From<LetoError>`).
//! - `first_derivative_simd` f32 SIMD extension that is CFD-specific.
//! - `differentiate_1d` / `differentiate_2d` / `laplacian_2d` free functions.

use super::schemes::FiniteDifferenceScheme;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

// ── Wrapper ──────────────────────────────────────────────────────────────────

/// Finite-difference operator for 1-D slice data.
///
/// Delegates to [`leto_ops::FiniteDifference`] (SSOT). Errors from the
/// underlying implementation are mapped through [`From<LetoError>`] to
/// [`cfd_core::Error`].
#[derive(Debug, Clone, Copy)]
pub struct FiniteDifference<T: RealField + Copy> {
    inner: leto_ops::FiniteDifference<T>,
}

impl<T: RealField + FloatElement + Copy> FiniteDifference<T> {
    /// Create with explicit scheme and grid spacing.
    #[must_use]
    pub fn new(scheme: FiniteDifferenceScheme, spacing: T) -> Self {
        Self {
            inner: leto_ops::FiniteDifference::new(scheme, spacing),
        }
    }

    /// Central-difference operator (2nd-order, default).
    #[must_use]
    pub fn central(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Central, spacing)
    }

    /// Forward-difference operator (1st-order).
    #[must_use]
    pub fn forward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Forward, spacing)
    }

    /// Backward-difference operator (1st-order).
    #[must_use]
    pub fn backward(spacing: T) -> Self {
        Self::new(FiniteDifferenceScheme::Backward, spacing)
    }

    /// Active scheme.
    #[must_use]
    pub fn scheme(&self) -> FiniteDifferenceScheme {
        self.inner.scheme()
    }

    /// Grid spacing.
    #[must_use]
    pub fn spacing(&self) -> T {
        self.inner.spacing()
    }

    /// Compute first derivative `f'` of `values`.
    ///
    /// # Errors
    /// - [`Error`] when fewer than 2 points are supplied.
    /// - [`Error`] when 2nd-order schemes need fewer than 3 points.
    pub fn first_derivative(&self, values: &[T]) -> Result<Array1<T>> {
        self.inner.first_derivative(values).map_err(Error::from)
    }

    /// Compute second derivative `f''` using central differences.
    ///
    /// # Errors
    /// Returns an error if the input array has fewer than 3 points.
    pub fn second_derivative(&self, values: &[T]) -> Result<Array1<T>> {
        self.inner.second_derivative(values).map_err(Error::from)
    }
}

impl<T: RealField + FloatElement + Copy + NumericElement> Default for FiniteDifference<T> {
    fn default() -> Self {
        Self::central(<T as NumericElement>::ONE)
    }
}

// ── f32 SIMD extension (CFD-specific) ────────────────────────────────────────

impl FiniteDifference<f32> {
    /// Compute first derivative using the f32 SIMD-friendly path.
    ///
    /// For `Central` scheme the inner loop uses hand-unrolled arithmetic that
    /// the compiler can autovectorize.  Other schemes fall back to the generic
    /// path.
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
        let inv_spacing = 1.0f32 / self.spacing();

        if self.scheme() == FiniteDifferenceScheme::Central {
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

// ── Free functions ────────────────────────────────────────────────────────────

/// Compute 1-D derivative using central differences.
///
/// # Errors
/// Returns an error if `values` has fewer than 2 points.
pub fn differentiate_1d<T: RealField + FloatElement + Copy>(
    values: &[T],
    spacing: T,
) -> Result<Array1<T>> {
    FiniteDifference::central(spacing).first_derivative(values)
}

/// Compute 2-D gradient (∂f/∂x, ∂f/∂y) using central differences.
///
/// # Errors
/// Returns an error if field dimensions are invalid.
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

/// Compute 2-D Laplacian ∇²f using central differences on a flat row-major slice.
///
/// # Errors
/// Returns an error if `field.len() != nx * ny`.
#[allow(clippy::similar_names)]
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
    let two = <T as FloatElement>::from_f64(2.0);

    for j in 1..ny - 1 {
        for i in 1..nx - 1 {
            let idx = j * nx + i;
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;
            let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Top and bottom boundaries
    for i in 0..nx {
        let idx = i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;
            let d2fdy2 = (field[idx + 2 * nx] - two * field[idx + nx] + field[idx]) / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }
        let idx = (ny - 1) * nx + i;
        if i > 0 && i < nx - 1 {
            let d2fdx2 = (field[idx + 1] - two * field[idx] + field[idx - 1]) / dx2;
            let d2fdy2 = (field[idx] - two * field[idx - nx] + field[idx - 2 * nx]) / dy2;
            laplacian[idx] = d2fdx2 + d2fdy2;
        }
    }

    // Left and right boundaries
    for j in 1..ny - 1 {
        let idx = j * nx;
        let d2fdx2 = (field[idx + 2] - two * field[idx + 1] + field[idx]) / dx2;
        let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;

        let idx = j * nx + (nx - 1);
        let d2fdx2 = (field[idx] - two * field[idx - 1] + field[idx - 2]) / dx2;
        let d2fdy2 = (field[idx + nx] - two * field[idx] + field[idx - nx]) / dy2;
        laplacian[idx] = d2fdx2 + d2fdy2;
    }

    Ok(laplacian)
}