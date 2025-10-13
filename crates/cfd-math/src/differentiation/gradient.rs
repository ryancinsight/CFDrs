//! Gradient computation for multi-dimensional fields.

use cfd_core::error::{Error, Result};
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Gradient computation for multi-dimensional fields
pub struct Gradient<T: RealField + Copy> {
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + From<f64> + FromPrimitive + Copy> Gradient<T> {
    /// Create gradient operator
    pub fn new(dx: T, dy: T, dz: T) -> Self {
        Self { dx, dy, dz }
    }

    /// Create uniform spacing gradient operator
    pub fn uniform(spacing: T) -> Self {
        Self::new(spacing, spacing, spacing)
    }

    /// Compute gradient of a 1D field
    /// 
    /// # Errors
    /// Returns an error if differentiation fails or field has insufficient points
    pub fn gradient_1d(&self, field: &[T]) -> Result<Vec<T>> {
        use super::FiniteDifference;
        let fd = FiniteDifference::central(self.dx);
        let grad = fd.first_derivative(field)?;
        // Zero-copy: convert DVector to Vec via into_iter (avoids clone)
        Ok(grad.iter().copied().collect())
    }

    /// Compute gradient of a 2D field (stored row-major)
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
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Use iterator combinators instead of nested loops
        gradients.extend((0..ny).flat_map(|j| {
            let dx = self.dx;
            let dy = self.dy;
            (0..nx).map(move |i| {
                let idx = j * nx + i;

                // Compute x-derivative
                let dfdx = if i == 0 {
                    // Forward difference
                    (field[idx + 1] - field[idx]) / dx
                } else if i == nx - 1 {
                    // Backward difference
                    (field[idx] - field[idx - 1]) / dx
                } else {
                    // Central difference
                    (field[idx + 1] - field[idx - 1]) / (two * dx)
                };

                // Compute y-derivative
                let dfdy = if j == 0 {
                    // Forward difference
                    (field[idx + nx] - field[idx]) / dy
                } else if j == ny - 1 {
                    // Backward difference
                    (field[idx] - field[idx - nx]) / dy
                } else {
                    // Central difference
                    (field[idx + nx] - field[idx - nx]) / (two * dy)
                };

                Vector3::new(dfdx, dfdy, T::zero())
            })
        }));

        Ok(gradients)
    }

    /// Compute gradient of a 3D field (stored in z-y-x order)
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
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        // Use iterator combinators for better performance
        gradients.extend((0..nz).flat_map(|k| {
            let dx = self.dx;
            let dy = self.dy;
            let dz = self.dz;
            (0..ny).flat_map(move |j| {
                let dx = dx;
                let dy = dy;
                let dz = dz;
                (0..nx).map(move |i| {
                    let idx = k * nx * ny + j * nx + i;

                    // Compute x-derivative
                    let dfdx = if i == 0 {
                        (field[idx + 1] - field[idx]) / dx
                    } else if i == nx - 1 {
                        (field[idx] - field[idx - 1]) / dx
                    } else {
                        (field[idx + 1] - field[idx - 1]) / (two * dx)
                    };

                    // Compute y-derivative
                    let dfdy = if j == 0 {
                        (field[idx + nx] - field[idx]) / dy
                    } else if j == ny - 1 {
                        (field[idx] - field[idx - nx]) / dy
                    } else {
                        (field[idx + nx] - field[idx - nx]) / (two * dy)
                    };

                    // Compute z-derivative
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

        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        let divergence: Vec<T> = (0..ny)
            .flat_map(|j| (0..nx).map(move |i| (i, j)))
            .map(|(i, j)| {
                let idx = j * nx + i;

                // Compute ∂u/∂x
                let dudx = if i == 0 {
                    (field[idx + 1].x - field[idx].x) / self.dx
                } else if i == nx - 1 {
                    (field[idx].x - field[idx - 1].x) / self.dx
                } else {
                    (field[idx + 1].x - field[idx - 1].x) / (two * self.dx)
                };

                // Compute ∂v/∂y
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
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Compute ∂v/∂x
                let dvdx = if i == 0 {
                    (field[idx + 1].y - field[idx].y) / self.dx
                } else if i == nx - 1 {
                    (field[idx].y - field[idx - 1].y) / self.dx
                } else {
                    (field[idx + 1].y - field[idx - 1].y) / (two * self.dx)
                };

                // Compute ∂u/∂y
                let dudy = if j == 0 {
                    (field[idx + nx].x - field[idx].x) / self.dy
                } else if j == ny - 1 {
                    (field[idx].x - field[idx - nx].x) / self.dy
                } else {
                    (field[idx + nx].x - field[idx - nx].x) / (two * self.dy)
                };

                // Curl in 2D: ∂v/∂x - ∂u/∂y
                curl.push(dvdx - dudy);
            }
        }

        Ok(curl)
    }
}

/// Compute gradient of a 2D field using central differences
pub fn compute_gradient_2d<T: RealField + From<f64> + FromPrimitive + Copy>(
    field: &[T],
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
) -> Result<Vec<Vector3<T>>> {
    Gradient::new(dx, dy, T::one()).gradient_2d(field, nx, ny)
}

/// Compute gradient of a 3D field using central differences
pub fn compute_gradient_3d<T: RealField + From<f64> + FromPrimitive + Copy>(
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
