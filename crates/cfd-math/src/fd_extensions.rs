//! CFD-specific finite-difference extensions.
//!
//! This module holds the parts of the old `differentiation` module that are
//! genuinely CFD-specific or higher-dimensional helpers built on top of the
//! leto-ops SSOT finite-difference surface:
//!
//! - `Gradient<T>` and its 2-D/3-D gradient helpers.
//! - `first_derivative_simd` f32 fast path.
//! - `differentiate_1d`, `differentiate_2d`, and `laplacian_2d` convenience
//!   functions.
//!
//! All generic 1-D finite-difference logic lives in `leto_ops::FiniteDifference`
//! (re-exported from `cfd_math::fd`); this module only provides the
//! CFD-layer convenience and vector-field helpers.

use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::geometry::Vector3;
use leto::Array1;
use leto_ops::{FiniteDifference, FiniteDifferenceScheme};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Compute the first derivative of a 1-D slice using central differences.
///
/// # Errors
/// Returns an error if fewer than 2 points are supplied.
pub fn differentiate_1d<T: RealField + FloatElement + Copy>(
    values: &[T],
    spacing: T,
) -> Result<Array1<T>> {
    FiniteDifference::central(spacing)
        .first_derivative(values)
        .map_err(Error::from)
}

/// Compute the -D gradient of a row-major field returning separate x/y components.
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
    let grad = Gradient::new(dx, dy, <T as NumericElement>::ONE);
    let gradients = grad.gradient_2d(field, nx, ny)?;

    let grad_x: Vec<T> = gradients.iter().map(|g| g.x).collect();
    let grad_y: Vec<T> = gradients.iter().map(|g| g.y).collect();

    Ok((grad_x, grad_y))
}

/// Compute the 2-D Laplacian ∇²f using central differences on a flat row-major slice.
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

/// Compute the first derivative using the f32 SIMD-friendly path.
///
/// For `Central` scheme the inner loop uses hand-unrolled arithmetic that
/// the compiler can autovectorize. Other schemes fall back to the generic
/// leto-ops path.
///
/// # Errors
/// Returns an error if the input array has fewer than 2 points.
pub fn first_derivative_simd(values: &[f32], spacing: f32, scheme: FiniteDifferenceScheme) -> Result<Vec<f32>> {
    if values.len() < 2 {
        return Err(Error::InvalidConfiguration(
            "Need at least 2 points for differentiation".to_string(),
        ));
    }

    let n = values.len();
    let mut result = vec![0.0f32; n];
    let inv_spacing = 1.0f32 / spacing;

    if scheme != FiniteDifferenceScheme::Central {
        let fd = FiniteDifference::new(scheme, spacing);
        return fd.first_derivative(values).map_err(Error::from).map(|arr| arr.iter().copied().collect());
    }

    if n > 2 {
        let scale = inv_spacing * 0.5;
        for i in 1..n - 1 {
            result[i] = (values[i + 1] - values[i - 1]) * scale;
        }
    }
    result[0] = (values[1] - values[0]) * inv_spacing;
    result[n - 1] = (values[n - 1] - values[n - 2]) * inv_spacing;

    Ok(result)
}

/// Gradient computation for multi-dimensional fields.
pub struct Gradient<T: RealField + Copy> {
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + FloatElement + Copy> Gradient<T> {
    /// Create gradient operator
    pub fn new(dx: T, dy: T, dz: T) -> Self {
        Self { dx, dy, dz }
    }

    /// Create uniform spacing gradient operator
    pub fn uniform(spacing: T) -> Self {
        Self::new(spacing, spacing, spacing)
    }

    /// Compute gradient of a 1-D field
    ///
    /// # Errors
    /// Returns an error if differentiation fails or field has insufficient points
    pub fn gradient_1d(&self, field: &[T]) -> Result<Vec<T>> {
        let fd = FiniteDifference::central(self.dx);
        let grad = fd.first_derivative(field).map_err(Error::from)?;
        // Copy the Leto result into the historical Vec return surface.
        Ok(grad.iter().copied().collect())
    }

    /// Compute gradient of a 2-D field (stored row-major)
    ///
    /// # Errors
    /// Returns an error if field dimensions don't match grid size or differentiation fails
    pub fn gradient_2d(&self, field: &[T], nx: usize, ny: usize) -> Result<Vec<Vector3<T>>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string(),
            ));
        }

        let mut gradients = Vec::with_capacity(nx * ny);
        let two = from_f64::<T>(2.0);

        gradients.extend((0..ny).flat_map(|j| {
            let dx = self.dx;
            let dy = self.dy;
            (0..nx).map(move |i| {
                let idx = j * nx + i;

                let dfdx = if i == 0 {
                    (field[idx + 1] - field[idx]) / dx
                } else if i == nx - 1 {
                    (field[idx] - field[idx - 1]) / dx
                } else {
                    (field[idx + 1] - field[idx - 1]) / (two * dx)
                };

                let dfdy = if j == 0 {
                    (field[idx + nx] - field[idx]) / dy
                } else if j == ny - 1 {
                    (field[idx] - field[idx - nx]) / dy
                } else {
                    (field[idx + nx] - field[idx - nx]) / (two * dy)
                };

                Vector3::new(dfdx, dfdy, <T as NumericElement>::ZERO)
            })
        }));

        Ok(gradients)
    }

    /// Compute gradient of a 3-D field (stored in z-y-x order)
    pub fn gradient_3d(
        &self,
        field: &[T],
        nx: usize,
        ny: usize,
        nz: usize,
    ) -> Result<Vec<Vector3<T>>> {
        if field.len() != nx * ny * nz {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string(),
            ));
        }

        let mut gradients = Vec::with_capacity(nx * ny * nz);
        let two = from_f64::<T>(2.0);

        gradients.extend((0..nz).flat_map(|k| {
            let dx = self.dx;
            let dy = self.dy;
            let dz = self.dz;
            (0..ny).flat_map(move |j| {
                (0..nx).map(move |i| {
                    let idx = k * nx * ny + j * nx + i;

                    let dfdx = if i == 0 {
                        (field[idx + 1] - field[idx]) / dx
                    } else if i == nx - 1 {
                        (field[idx] - field[idx - 1]) / dx
                    } else {
                        (field[idx + 1] - field[idx - 1]) / (two * dx)
                    };

                    let dfdy = if j == 0 {
                        (field[idx + nx] - field[idx]) / dy
                    } else if j == ny - 1 {
                        (field[idx] - field[idx - nx]) / dy
                    } else {
                        (field[idx + nx] - field[idx - nx]) / (two * dy)
                    };

                    let dfdz = if k == 0 {
                        (field[idx + nx * ny] - field[idx]) / dz
                    } else if k == nz - 1 {
                        (field[idx] - field[idx - nx * ny]) / dz
                    } else {
                        (field[idx + nx * ny] - field[idx - nx * ny]) / (two * dz)
                    };

                    Vector3::new(dfdx, dfdy, dfdz)
                })
            })
        }));

        Ok(gradients)
    }

    /// Compute divergence of a vector field in 2D
    pub fn divergence_2d(&self, field: &[Vector3<T>], nx: usize, ny: usize) -> Result<Vec<T>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string(),
            ));
        }

        let two = from_f64::<T>(2.0);

        let divergence: Vec<T> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| {
                let idx = j * nx + i;

                let dudx = if i == 0 {
                    (field[idx + 1].x - field[idx].x) / self.dx
                } else if i == nx - 1 {
                    (field[idx].x - field[idx - 1].x) / self.dx
                } else {
                    (field[idx + 1].x - field[idx - 1].x) / (two * self.dx)
                };

                let dvdy = if j == 0 {
                    (field[idx + nx].y - field[idx].y) / self.dy
                } else if j == ny - 1 {
                    (field[idx].y - field[idx - nx].y) / self.dy
                } else {
                    (field[idx + nx].y - field[idx - nx].y) / (two * self.dy)
                };

                dudx + dvdy
            })
            .collect();

        Ok(divergence)
    }

    /// Compute curl of a vector field in 2D (returns z-component only)
    pub fn curl_2d(&self, field: &[Vector3<T>], nx: usize, ny: usize) -> Result<Vec<T>> {
        if field.len() != nx * ny {
            return Err(Error::InvalidConfiguration(
                "Field size doesn't match grid dimensions".to_string(),
            ));
        }

        let mut curl = Vec::with_capacity(nx * ny);
        let two = from_f64::<T>(2.0);

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                let dvdx = if i == 0 {
                    (field[idx + 1].y - field[idx].y) / self.dx
                } else if i == nx - 1 {
                    (field[idx].y - field[idx - 1].y) / self.dx
                } else {
                    (field[idx + 1].y - field[idx - 1].y) / (two * self.dx)
                };

                let dudy = if j == 0 {
                    (field[idx + nx].x - field[idx].x) / self.dy
                } else if j == ny - 1 {
                    (field[idx].x - field[idx - nx].x) / self.dy
                } else {
                    (field[idx + nx].x - field[idx - nx].x) / (two * self.dy)
                };

                curl.push(dvdx - dudy);
            }
        }

        Ok(curl)
    }
}

/// Compute gradient of a 2D field using central differences
pub fn compute_gradient_2d<T: RealField + FloatElement + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Result<Vec<Vector3<T>>> {
    Gradient::new(dx, dy, <T as NumericElement>::ONE).gradient_2d(field, nx, ny)
}

/// Compute gradient of a 3D field using central differences
pub fn compute_gradient_3d<T: RealField + FloatElement + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
) -> Result<Vec<Vector3<T>>> {
    Gradient::new(dx, dy, dz).gradient_3d(field, nx, ny, nz)
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_gradient_2d() {
        let grad = Gradient::uniform(1.0);

        let nx = 3;
        let ny = 3;
        let mut field = vec![0.0; nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = x * x + y * y;
            }
        }

        let gradients = grad
            .gradient_2d(&field, nx, ny)
            .expect("Failed to compute gradient in test");

        let center_grad = &gradients[nx + 1];
        assert_relative_eq!(center_grad.x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.y, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_divergence_2d() {
        let grad = Gradient::uniform(1.0);

        let nx = 3;
        let ny = 3;
        let mut field = vec![Vector3::zeros(); nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = Vector3::new(x, y, 0.0);
            }
        }

        let divergence = grad
            .divergence_2d(&field, nx, ny)
            .expect("Failed to compute divergence in test");

        for &div in &divergence {
            assert_relative_eq!(div, 2.0, epsilon = 1e-10);
        }
    }

    #[test]
    fn test_curl_2d() {
        let grad = Gradient::uniform(1.0);

        let nx = 3;
        let ny = 3;
        let mut field = vec![Vector3::zeros(); nx * ny];

        for j in 0..ny {
            for i in 0..nx {
                let x = i as f64;
                let y = j as f64;
                field[j * nx + i] = Vector3::new(-y, x, 0.0);
            }
        }

        let curl = grad
            .curl_2d(&field, nx, ny)
            .expect("Failed to compute curl in test");

        assert_relative_eq!(curl[nx + 1], 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_gradient_3d() {
        let grad = Gradient::uniform(1.0);

        let nx = 3;
        let ny = 3;
        let nz = 3;
        let mut field = vec![0.0; nx * ny * nz];

        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    let x = i as f64;
                    let y = j as f64;
                    let z = k as f64;
                    field[k * nx * ny + j * nx + i] = x * x + y * y + z * z;
                }
            }
        }

        let gradients = grad
            .gradient_3d(&field, nx, ny, nz)
            .expect("Failed to compute gradient in test");

        let center_grad = &gradients[nx * ny + nx + 1];
        assert_relative_eq!(center_grad.x, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.y, 2.0, epsilon = 1e-10);
        assert_relative_eq!(center_grad.z, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_first_derivative_simd_central() {
        let values: Vec<f32> = (0..=10).map(|i| (i as f32).powi(2)).collect();
        let derivatives = first_derivative_simd(&values, 1.0, FiniteDifferenceScheme::Central)
            .expect("SIMD derivative should succeed");

        // For f(x)=x^2, central derivative at interior points is ~2x.
        assert_relative_eq!(derivatives[2], 4.0, epsilon = 1e-5);
        assert_relative_eq!(derivatives[5], 10.0, epsilon = 1e-5);
    }
}
