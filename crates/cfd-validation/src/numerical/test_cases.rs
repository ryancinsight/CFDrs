//! Test case generation for numerical validation

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::FromPrimitive;

// Named constants for test cases
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/// Create diagonal system for testing
pub fn create_diagonal_system<T: RealField + Copy + FromPrimitive>(
    n: usize,
) -> Result<(CsrMatrix<T>, DVector<T>, DVector<T>)> {
    // Create diagonal matrix with entries 1, 2, 3, ..., n using iterators
    let diagonal_entries: Vec<(usize, usize, T)> = (0..n)
        .map(|i| (i, i, T::from_usize(i + 1).unwrap_or_else(T::zero)))
        .collect();

    let (row_indices, col_indices, values): (Vec<_>, Vec<_>, Vec<_>) =
        diagonal_entries.into_iter().fold(
            (Vec::new(), Vec::new(), Vec::new()),
            |(mut rows, mut cols, mut vals), (r, c, v)| {
                rows.push(r);
                cols.push(c);
                vals.push(v);
                (rows, cols, vals)
            },
        );

    let a = build_csr_matrix(n, n, row_indices, col_indices, values)?;

    // Create RHS b = [1, 1, ..., 1]
    let b = DVector::from_element(n, T::one());

    // Analytical solution: x[i] = 1 / (i + 1)
    let analytical = DVector::from_iterator(
        n,
        (0..n).map(|i| T::one() / T::from_usize(i + 1).unwrap_or_else(T::zero)),
    );

    Ok((a, b, analytical))
}

/// Create tridiagonal system (1D Poisson)
pub fn create_tridiagonal_system<T: RealField + Copy + FromPrimitive>(
    n: usize,
) -> Result<(CsrMatrix<T>, DVector<T>, DVector<T>)> {
    let h = T::one() / T::from_usize(n + 1).unwrap_or_else(T::zero);
    let h_squared = h * h;

    // Create tridiagonal matrix for -u'' using iterators
    let diagonal_value = T::from_f64(TWO).unwrap_or_else(T::zero) / h_squared;
    let off_diagonal_value = -T::one() / h_squared;

    let (row_indices, col_indices, values): (Vec<_>, Vec<_>, Vec<_>) = (0..n)
        .flat_map(|i| {
            let mut entries = vec![(i, i, diagonal_value)];

            if i > 0 {
                entries.push((i, i - 1, off_diagonal_value));
            }
            if i < n - 1 {
                entries.push((i, i + 1, off_diagonal_value));
            }

            entries
        })
        .fold(
            (Vec::new(), Vec::new(), Vec::new()),
            |(mut rows, mut cols, mut vals), (r, c, v)| {
                rows.push(r);
                cols.push(c);
                vals.push(v);
                (rows, cols, vals)
            },
        );

    let a = build_csr_matrix(n, n, row_indices, col_indices, values)?;

    // RHS for manufactured solution u(x) = x(1-x)
    let b = DVector::from_iterator(
        n,
        (1..=n).map(|_| {
            // For the manufactured solution u(x) = x(1-x), the Laplacian is -2
            T::from_f64(TWO).unwrap_or_else(T::zero)
        }),
    );

    // Analytical solution
    let analytical = DVector::from_iterator(
        n,
        (1..=n).map(|i| {
            let x = T::from_usize(i).unwrap_or_else(T::zero) * h;
            x * (T::one() - x)
        }),
    );

    Ok((a, b, analytical))
}

/// Create 2D Poisson system with 5-point stencil discretization
/// Solves: -∇²u = f on unit square with Dirichlet boundary conditions
pub fn create_2d_poisson_system<T: RealField + Copy + FromPrimitive>(
    nx: usize,
    ny: usize,
) -> Result<(CsrMatrix<T>, DVector<T>, DVector<T>)> {
    let n = nx * ny;
    let h = T::one() / T::from_usize(nx - 1).unwrap_or_else(T::zero);
    let h2 = h * h;

    // Build sparse matrix using 5-point stencil
    let mut row_offsets = vec![0];
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    // Process each grid point
    for idx in 0..n {
        let i = idx % nx;
        let j = idx / nx;

        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            // Boundary point: u = 0 (identity row)
            col_indices.push(idx);
            values.push(T::one());
        } else {
            // Interior point: 5-point stencil
            let center_coeff = T::from_f64(FOUR).unwrap_or_else(T::zero) / h2;
            let neighbor_coeff = -T::one() / h2;

            // Left neighbor
            col_indices.push(idx - 1);
            values.push(neighbor_coeff);

            // Bottom neighbor
            col_indices.push(idx - nx);
            values.push(neighbor_coeff);

            // Center
            col_indices.push(idx);
            values.push(center_coeff);

            // Top neighbor
            col_indices.push(idx + nx);
            values.push(neighbor_coeff);

            // Right neighbor
            col_indices.push(idx + 1);
            values.push(neighbor_coeff);
        }

        row_offsets.push(col_indices.len());
    }

    let a = CsrMatrix::try_from_csr_data(n, n, row_offsets, col_indices, values)
        .map_err(|_| Error::InvalidConfiguration("Matrix construction failed".into()))?;

    // Create RHS with manufactured solution u(x,y) = sin(πx)sin(πy)
    let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::zero);
    let mut b = DVector::zeros(n);
    let mut analytical = DVector::zeros(n);

    for idx in 0..n {
        let i = idx % nx;
        let j = idx / nx;
        let x = T::from_usize(i).unwrap_or_else(T::zero) * h;
        let y = T::from_usize(j).unwrap_or_else(T::zero) * h;

        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            b[idx] = T::zero();
            analytical[idx] = T::zero();
        } else {
            // f = 2π²sin(πx)sin(πy)
            b[idx] = T::from_f64(TWO).unwrap_or_else(T::zero)
                * pi
                * pi
                * (pi * x).sin()
                * (pi * y).sin();
            analytical[idx] = (pi * x).sin() * (pi * y).sin();
        }
    }

    Ok((a, b, analytical))
}

/// Create Hilbert matrix system (ill-conditioned)
pub fn create_hilbert_system<T: RealField + Copy + FromPrimitive>(
    n: usize,
) -> Result<(CsrMatrix<T>, DVector<T>, DVector<T>)> {
    let mut row_indices = Vec::new();
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    // Create Hilbert matrix H[i,j] = 1/(i+j+1)
    for i in 0..n {
        for j in 0..n {
            row_indices.push(i);
            col_indices.push(j);
            values.push(T::one() / T::from_usize(i + j + 1).unwrap_or_else(T::zero));
        }
    }

    let a = build_csr_matrix(n, n, row_indices, col_indices, values)?;

    // Create RHS with known solution [1, 1, ..., 1]
    let analytical = DVector::from_element(n, T::one());

    // Compute b = A * x_true
    let b = &a * &analytical;

    Ok((a, b, analytical))
}

/// Helper function to build CSR matrix from triplets
fn build_csr_matrix<T: RealField + Copy>(
    rows: usize,
    cols: usize,
    row_indices: Vec<usize>,
    col_indices: Vec<usize>,
    values: Vec<T>,
) -> Result<CsrMatrix<T>> {
    // Convert triplets to CSR format
    let mut row_offsets = vec![0; rows + 1];
    let mut sorted_data: Vec<(usize, usize, T)> = row_indices
        .into_iter()
        .zip(col_indices)
        .zip(values)
        .map(|((r, c), v)| (r, c, v))
        .collect();

    sorted_data.sort_by_key(|(r, c, _)| (*r, *c));

    let mut current_row = 0;
    let mut csr_col_indices = Vec::new();
    let mut csr_values = Vec::new();

    for (r, c, v) in sorted_data {
        while current_row < r {
            current_row += 1;
            row_offsets[current_row] = csr_col_indices.len();
        }
        csr_col_indices.push(c);
        csr_values.push(v);
    }

    while current_row < rows {
        current_row += 1;
        row_offsets[current_row] = csr_col_indices.len();
    }

    CsrMatrix::try_from_csr_data(rows, cols, row_offsets, csr_col_indices, csr_values)
        .map_err(|_| Error::InvalidConfiguration("Matrix construction failed".into()))
}
