//! Test case generation for numerical validation

use crate::scalar;
use cfd_core::error::{Error, Result};
use cfd_math::sparse::SparseMatrix as CsrMatrix;
use eunomia::{FloatElement, RealField};
use leto::Array1;
use leto_ops::Scalar as LetoScalar;

// Named constants for test cases
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/// Create diagonal system for testing
pub fn create_diagonal_system<T: RealField + Copy + FloatElement + LetoScalar>(
    n: usize,
) -> Result<(CsrMatrix<T>, Array1<T>, Array1<T>)> {
    // Create diagonal matrix with entries 1, 2, 3, ..., n using iterators
    let diagonal_entries: Vec<(usize, usize, T)> =
        (0..n).map(|i| (i, i, scalar::from_usize(i + 1))).collect();

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
    let b = Array1::from_elem([n], scalar::one::<T>());

    // Analytical solution: x[i] = 1 / (i + 1)
    let analytical = Array1::from_shape_vec(
        [n],
        (0..n)
            .map(|i| scalar::one::<T>() / scalar::from_usize(i + 1))
            .collect(),
    )
    .map_err(|error| Error::InvalidConfiguration(error.to_string()))?;

    Ok((a, b, analytical))
}

/// Create tridiagonal system (1D Poisson)
pub fn create_tridiagonal_system<T: RealField + Copy + FloatElement + LetoScalar>(
    n: usize,
) -> Result<(CsrMatrix<T>, Array1<T>, Array1<T>)> {
    let h = scalar::one::<T>() / scalar::from_usize::<T>(n + 1);
    let h_squared = h * h;

    // Create tridiagonal matrix for -u'' using iterators
    let diagonal_value = scalar::from_f64::<T>(TWO) / h_squared;
    let off_diagonal_value = -scalar::one::<T>() / h_squared;

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
            (Vec::<usize>::new(), Vec::<usize>::new(), Vec::<T>::new()),
            |(mut rows, mut cols, mut vals), (r, c, v)| {
                rows.push(r);
                cols.push(c);
                vals.push(v);
                (rows, cols, vals)
            },
        );

    let a = build_csr_matrix(n, n, row_indices, col_indices, values)?;

    // RHS for manufactured solution u(x) = x(1-x)
    let b = Array1::from_shape_vec(
        [n],
        (1..=n)
            .map(|_| {
                // For the manufactured solution u(x) = x(1-x), the Laplacian is -2
                scalar::from_f64(TWO)
            })
            .collect(),
    );
    let b = b.map_err(|error| Error::InvalidConfiguration(error.to_string()))?;

    // Analytical solution
    let analytical = Array1::from_shape_vec(
        [n],
        (1..=n)
            .map(|i| {
                let x = scalar::from_usize::<T>(i) * h;
                x * (scalar::one::<T>() - x)
            })
            .collect(),
    )
    .map_err(|error| Error::InvalidConfiguration(error.to_string()))?;

    Ok((a, b, analytical))
}

/// Create 2D Poisson system with 5-point stencil discretization
/// Solves: -∇²u = f on unit square with Dirichlet boundary conditions
pub fn create_2d_poisson_system<T: RealField + Copy + FloatElement + LetoScalar>(
    nx: usize,
    ny: usize,
) -> Result<(CsrMatrix<T>, Array1<T>, Array1<T>)> {
    let n = nx * ny;
    let h = scalar::one::<T>() / scalar::from_usize::<T>(nx - 1);
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
            values.push(scalar::one::<T>());
        } else {
            // Interior point: 5-point stencil
            let center_coeff = scalar::from_f64::<T>(FOUR) / h2;
            let neighbor_coeff = -scalar::one::<T>() / h2;

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

    let a = CsrMatrix::from_parts(values, col_indices, row_offsets, n, n).map_err(|error| {
        Error::InvalidConfiguration(format!("Matrix construction failed: {error}"))
    })?;

    // Create RHS with manufactured solution u(x,y) = sin(πx)sin(πy)
    let pi = scalar::from_f64::<T>(std::f64::consts::PI);
    let mut b = Array1::zeros([n]);
    let mut analytical = Array1::zeros([n]);

    for idx in 0..n {
        let i = idx % nx;
        let j = idx / nx;
        let x = scalar::from_usize::<T>(i) * h;
        let y = scalar::from_usize::<T>(j) * h;

        if i == 0 || i == nx - 1 || j == 0 || j == ny - 1 {
            b[idx] = scalar::zero::<T>();
            analytical[idx] = scalar::zero::<T>();
        } else {
            // f = 2π²sin(πx)sin(πy)
            let sin_x = scalar::sin(pi * x);
            let sin_y = scalar::sin(pi * y);
            b[idx] = scalar::from_f64::<T>(TWO) * pi * pi * sin_x * sin_y;
            analytical[idx] = sin_x * sin_y;
        }
    }

    Ok((a, b, analytical))
}

/// Create Hilbert matrix system (ill-conditioned)
pub fn create_hilbert_system<T: RealField + Copy + FloatElement + LetoScalar>(
    n: usize,
) -> Result<(CsrMatrix<T>, Array1<T>, Array1<T>)> {
    let mut row_indices = Vec::new();
    let mut col_indices = Vec::new();
    let mut values = Vec::new();

    // Create Hilbert matrix H[i,j] = 1/(i+j+1)
    for i in 0..n {
        for j in 0..n {
            row_indices.push(i);
            col_indices.push(j);
            values.push(scalar::one::<T>() / scalar::from_usize::<T>(i + j + 1));
        }
    }

    let a = build_csr_matrix(n, n, row_indices, col_indices, values)?;

    // Create RHS with known solution [1, 1, ..., 1]
    let analytical = Array1::from_elem([n], scalar::one::<T>());

    // Compute b = A * x_true
    let mut b = Array1::zeros([n]);
    cfd_math::sparse::try_spmv(&a, &analytical, &mut b)?;

    Ok((a, b, analytical))
}

/// Helper function to build CSR matrix from triplets
fn build_csr_matrix<T: RealField + Copy + LetoScalar>(
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

    CsrMatrix::from_parts(csr_values, csr_col_indices, row_offsets, rows, cols).map_err(|error| {
        Error::InvalidConfiguration(format!("Matrix construction failed: {error}"))
    })
}
