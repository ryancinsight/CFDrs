//! Leto dense-solve bridge for CSR storage boundaries.

use cfd_core::error::{Error, Result};
use leto::{Array1, Array2};
use leto_ops::{solve as leto_lu_solve, CsrMatrix as LetoCsrMatrix, RealScalar as LetoRealScalar};

pub(super) fn solve_leto_csr_with_leto_dense_array<T>(
    matrix: &LetoCsrMatrix<T>,
    rhs: &Array1<T>,
) -> Result<Array1<T>>
where
    T: LetoRealScalar,
{
    let nrows = matrix.nrows();
    let ncols = matrix.ncols();
    let mut dense_values = vec![T::ZERO; nrows * ncols];

    let row_offsets = matrix.row_ptr();
    let col_indices = matrix.col_indices();
    let values = matrix.values();

    for row in 0..nrows {
        for idx in row_offsets[row]..row_offsets[row + 1] {
            dense_values[row * ncols + col_indices[idx]] = values[idx];
        }
    }

    let dense = Array2::from_shape_vec([nrows, ncols], dense_values).map_err(|error| {
        Error::InvalidConfiguration(format!("Invalid Leto dense fallback matrix: {error}"))
    })?;

    leto_lu_solve(&dense.view(), &rhs.view())
        .map_err(|error| Error::Solver(format!("Dense LU failed through Leto provider: {error}")))
}
