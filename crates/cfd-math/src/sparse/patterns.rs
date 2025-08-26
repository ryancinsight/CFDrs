//! Common sparse matrix patterns for CFD

use super::builder::SparseMatrixBuilder;
use cfd_core::error::{Error, Result};
use nalgebra::RealField;
use nalgebra_sparse::CsrMatrix;

/// Factory for common sparse matrix patterns
pub struct SparsePatterns;

impl SparsePatterns {
    /// Create tridiagonal matrix pattern
    pub fn tridiagonal<T>(n: usize, lower: T, diag: T, upper: T) -> Result<CsrMatrix<T>>
    where
        T: RealField + Copy,
    {
        if n == 0 {
            return Err(Error::InvalidConfiguration(
                "Matrix size must be positive".to_string(),
            ));
        }

        let mut builder = SparseMatrixBuilder::with_capacity(n, n, 3 * n - 2);

        // First row
        builder.add_entry(0, 0, diag)?;
        if n > 1 {
            builder.add_entry(0, 1, upper)?;
        }

        // Middle rows
        for i in 1..n - 1 {
            builder.add_entry(i, i - 1, lower)?;
            builder.add_entry(i, i, diag)?;
            builder.add_entry(i, i + 1, upper)?;
        }

        // Last row
        if n > 1 {
            builder.add_entry(n - 1, n - 2, lower)?;
            builder.add_entry(n - 1, n - 1, diag)?;
        }

        builder.build()
    }

    /// Create five-point stencil for 2D Laplacian
    pub fn five_point_stencil<T>(nx: usize, ny: usize, dx: T, dy: T) -> Result<CsrMatrix<T>>
    where
        T: RealField + Copy,
    {
        if nx == 0 || ny == 0 {
            return Err(Error::InvalidConfiguration(
                "Grid dimensions must be positive".to_string(),
            ));
        }

        let n = nx * ny;
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let cx = T::one() / dx2;
        let cy = T::one() / dy2;
        let center = -(T::from_f64(2.0).unwrap_or(T::one() + T::one()) * (cx + cy));

        let mut builder = SparseMatrixBuilder::with_capacity(n, n, 5 * n);

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Center point
                builder.add_entry(idx, idx, center)?;

                // West neighbor
                if i > 0 {
                    builder.add_entry(idx, idx - 1, cx)?;
                }

                // East neighbor
                if i < nx - 1 {
                    builder.add_entry(idx, idx + 1, cx)?;
                }

                // South neighbor
                if j > 0 {
                    builder.add_entry(idx, idx - nx, cy)?;
                }

                // North neighbor
                if j < ny - 1 {
                    builder.add_entry(idx, idx + nx, cy)?;
                }
            }
        }

        builder.build()
    }

    /// Create nine-point stencil for 2D biharmonic operator
    pub fn nine_point_stencil<T>(nx: usize, ny: usize, dx: T, dy: T) -> Result<CsrMatrix<T>>
    where
        T: RealField + Copy,
    {
        if nx < 3 || ny < 3 {
            return Err(Error::InvalidConfiguration(
                "Grid must be at least 3x3 for nine-point stencil".to_string(),
            ));
        }

        let n = nx * ny;
        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let dxdy = dx * dy;

        // Stencil coefficients
        let corner = T::one() / (T::from_f64(4.0).unwrap_or(T::one()) * dxdy);
        let x_side = T::one() / dx2;
        let y_side = T::one() / dy2;
        let center = -(T::from_f64(2.0).unwrap_or(T::one() + T::one()) * (x_side + y_side)
            + T::from_f64(4.0).unwrap_or(T::one()) * corner);

        let mut builder = SparseMatrixBuilder::with_capacity(n, n, 9 * n);

        for j in 0..ny {
            for i in 0..nx {
                let idx = j * nx + i;

                // Center
                builder.add_entry(idx, idx, center)?;

                // Cardinal directions
                if i > 0 {
                    builder.add_entry(idx, idx - 1, x_side)?;
                }
                if i < nx - 1 {
                    builder.add_entry(idx, idx + 1, x_side)?;
                }
                if j > 0 {
                    builder.add_entry(idx, idx - nx, y_side)?;
                }
                if j < ny - 1 {
                    builder.add_entry(idx, idx + nx, y_side)?;
                }

                // Diagonal directions
                if i > 0 && j > 0 {
                    builder.add_entry(idx, idx - nx - 1, corner)?;
                }
                if i < nx - 1 && j > 0 {
                    builder.add_entry(idx, idx - nx + 1, corner)?;
                }
                if i > 0 && j < ny - 1 {
                    builder.add_entry(idx, idx + nx - 1, corner)?;
                }
                if i < nx - 1 && j < ny - 1 {
                    builder.add_entry(idx, idx + nx + 1, corner)?;
                }
            }
        }

        builder.build()
    }

    /// Create banded matrix pattern
    pub fn banded<T>(n: usize, bandwidth: usize, value: T) -> Result<CsrMatrix<T>>
    where
        T: RealField + Copy,
    {
        if n == 0 || bandwidth == 0 {
            return Err(Error::InvalidConfiguration(
                "Matrix size and bandwidth must be positive".to_string(),
            ));
        }

        let mut builder = SparseMatrixBuilder::new(n, n);

        for i in 0..n {
            for j in i.saturating_sub(bandwidth)..=(i + bandwidth).min(n - 1) {
                builder.add_entry(i, j, value)?;
            }
        }

        builder.build()
    }
}
