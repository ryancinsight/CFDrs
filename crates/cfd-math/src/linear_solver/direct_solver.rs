//! Direct sparse solver for linear systems.
//!
//! Uses rsparse (CSparse-based LU) for robust direct solves of sparse matrices.

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use rsparse::data::Sprs;

/// Direct sparse solver configuration.
#[derive(Debug, Clone)]
pub struct DirectSparseSolver {
    /// Maximum system size (number of rows) to attempt with direct solve.
    pub max_size: usize,
    /// Ordering strategy: -1 natural, 0 Cholesky, 1 LU, 2 QR.
    pub ordering: i8,
    /// Pivot tolerance for LU factorization.
    pub pivot_tolerance: f64,
}

impl Default for DirectSparseSolver {
    fn default() -> Self {
        Self {
            max_size: 20000,
            ordering: 1,
            pivot_tolerance: 1e-12,
        }
    }
}

impl DirectSparseSolver {
    /// Returns true if the solver should attempt a direct solve for this size.
    pub fn can_handle_size(&self, size: usize) -> bool {
        size <= self.max_size
    }

    /// Solve Ax = b using sparse LU factorization.
    pub fn solve<T>(&self, matrix: &CsrMatrix<T>, rhs: &DVector<T>) -> Result<DVector<T>>
    where
        T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
    {
        let nrows = matrix.nrows();
        let ncols = matrix.ncols();

        if nrows != ncols {
            return Err(Error::InvalidConfiguration(
                "Direct solver requires a square matrix".to_string(),
            ));
        }
        if rhs.len() != nrows {
            return Err(Error::InvalidConfiguration(
                "RHS dimension does not match matrix size".to_string(),
            ));
        }

        let (sprs, mut b) = csr_to_csc_sprs_f64(matrix, rhs)?;

        rsparse::lusol(&sprs, &mut b, self.ordering, self.pivot_tolerance)
            .map_err(|e| Error::Solver(format!("Direct sparse LU failed: {}", e)))?;

        let mut x = DVector::zeros(nrows);
        for (i, value) in b.iter().enumerate() {
            x[i] = T::from_f64(*value).unwrap_or_else(T::zero);
        }
        Ok(x)
    }
}

fn csr_to_csc_sprs_f64<T>(matrix: &CsrMatrix<T>, rhs: &DVector<T>) -> Result<(Sprs<f64>, Vec<f64>)>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive,
{
    let nrows = matrix.nrows();
    let ncols = matrix.ncols();
    let nnz = matrix.nnz();
    let row_offsets = matrix.row_offsets();
    let col_indices = matrix.col_indices();
    let values = matrix.values();

    let mut col_counts = vec![0usize; ncols];
    for &col in col_indices {
        col_counts[col] += 1;
    }

    let mut p = vec![0isize; ncols + 1];
    for col in 0..ncols {
        p[col + 1] = p[col] + col_counts[col] as isize;
    }

    let mut next = p[..ncols].iter().map(|&v| v as usize).collect::<Vec<_>>();
    let mut i = vec![0usize; nnz];
    let mut x = vec![0f64; nnz];

    for row in 0..nrows {
        let row_start = row_offsets[row];
        let row_end = row_offsets[row + 1];
        for idx in row_start..row_end {
            let col = col_indices[idx];
            let dst = next[col];
            i[dst] = row;
            let value = values[idx].to_f64().ok_or_else(|| {
                Error::InvalidConfiguration("Matrix value not convertible to f64".to_string())
            })?;
            x[dst] = value;
            next[col] += 1;
        }
    }

    let mut b = Vec::with_capacity(rhs.len());
    for value in rhs.iter() {
        b.push(value.to_f64().ok_or_else(|| {
            Error::InvalidConfiguration("RHS value not convertible to f64".to_string())
        })?);
    }

    let sprs = Sprs {
        nzmax: nnz,
        m: nrows,
        n: ncols,
        x,
        i,
        p,
    };

    Ok((sprs, b))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sparse::SparseMatrixBuilder;

    #[test]
    fn direct_solver_solves_small_system() {
        let mut builder = SparseMatrixBuilder::new(2, 2);
        builder.add_entry(0, 0, 3.0).unwrap();
        builder.add_entry(0, 1, 1.0).unwrap();
        builder.add_entry(1, 0, 1.0).unwrap();
        builder.add_entry(1, 1, 2.0).unwrap();

        let mut rhs = DVector::from_vec(vec![9.0, 8.0]);
        let matrix = builder.build_with_rhs(&mut rhs).unwrap();

        let solver = DirectSparseSolver::default();
        let x = solver.solve(&matrix, &rhs).unwrap();

        assert!((x[0] - 2.0).abs() < 1e-8);
        assert!((x[1] - 3.0).abs() < 1e-8);
    }
}
