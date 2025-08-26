//! Sparse matrix operations and extensions

use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};
use nalgebra_sparse::CsrMatrix;
use num_traits::{Float, FromPrimitive, Signed};

/// Extension trait for sparse matrix operations
pub trait SparseMatrixExt<T: RealField + Copy> {
    /// Extract diagonal elements
    fn diagonal(&self) -> DVector<T>;

    /// Set diagonal elements
    fn set_diagonal(&mut self, diag: &DVector<T>) -> Result<()>;

    /// Scale matrix by a scalar
    fn scale(&mut self, factor: T);

    /// Add identity matrix scaled by factor
    fn add_identity(&mut self, factor: T) -> Result<()>;

    /// Compute Frobenius norm
    fn frobenius_norm(&self) -> T
    where
        T: Float;

    /// Compute condition number estimate
    fn condition_estimate(&self) -> Result<T>
    where
        T: Float + FromPrimitive;

    /// Check if matrix is diagonally dominant
    fn is_diagonally_dominant(&self) -> bool;

    /// Apply row scaling
    fn scale_rows(&mut self, scaling: &DVector<T>) -> Result<()>;

    /// Apply column scaling  
    fn scale_columns(&mut self, scaling: &DVector<T>) -> Result<()>;
}

impl<T: RealField + Copy> SparseMatrixExt<T> for CsrMatrix<T> {
    fn diagonal(&self) -> DVector<T> {
        let mut diag = DVector::zeros(self.nrows().min(self.ncols()));

        for i in 0..diag.len() {
            // Access row and find diagonal element
            let row = self.row(i);
            for (col_idx, &col) in row.col_indices().iter().enumerate() {
                if col == i {
                    diag[i] = row.values()[col_idx];
                    break;
                }
            }
        }

        diag
    }

    fn set_diagonal(&mut self, diag: &DVector<T>) -> Result<()> {
        if diag.len() != self.nrows().min(self.ncols()) {
            return Err(Error::InvalidConfiguration(
                "Diagonal size mismatch".to_string(),
            ));
        }

        // This requires mutable access to CSR internals
        // For now, we'll return an error as CSR matrices are not easily mutable
        Err(Error::InvalidConfiguration(
            "Direct diagonal modification not supported for CSR format".to_string(),
        ))
    }

    fn scale(&mut self, factor: T) {
        // Scale all values in the matrix
        for value in self.values_mut() {
            *value *= factor;
        }
    }

    fn add_identity(&mut self, _factor: T) -> Result<()> {
        if self.nrows() != self.ncols() {
            return Err(Error::InvalidConfiguration(
                "Matrix must be square to add identity".to_string(),
            ));
        }

        // This operation requires rebuilding the matrix structure
        // For CSR format, this is complex and would require conversion
        Err(Error::InvalidConfiguration(
            "Adding identity not directly supported for CSR format".to_string(),
        ))
    }

    fn frobenius_norm(&self) -> T
    where
        T: Float,
    {
        Float::sqrt(
            self.values()
                .iter()
                .map(|&v| v * v)
                .fold(T::zero(), |acc, v| acc + v),
        )
    }

    fn condition_estimate(&self) -> Result<T>
    where
        T: Float + FromPrimitive,
    {
        if self.nrows() != self.ncols() {
            return Err(Error::InvalidConfiguration(
                "Condition number requires square matrix".to_string(),
            ));
        }

        // Estimate using diagonal dominance
        let diag = self.diagonal();
        let mut max_ratio = T::one();

        for i in 0..self.nrows() {
            if Signed::abs(&diag[i]) < T::from_f64(1e-12).unwrap_or(T::epsilon()) {
                return Ok(T::infinity());
            }

            let row = self.row(i);
            let row_sum: T = row
                .values()
                .iter()
                .enumerate()
                .filter(|(idx, _)| row.col_indices()[*idx] != i)
                .map(|(_, &v)| Signed::abs(&v))
                .fold(T::zero(), |acc, v| acc + v);

            let ratio = (row_sum + Signed::abs(&diag[i])) / Signed::abs(&diag[i]);
            if ratio > max_ratio {
                max_ratio = ratio;
            }
        }

        Ok(max_ratio)
    }

    fn is_diagonally_dominant(&self) -> bool {
        if self.nrows() != self.ncols() {
            return false;
        }

        for i in 0..self.nrows() {
            let row = self.row(i);
            let mut diag_val = T::zero();
            let mut off_diag_sum = T::zero();

            for (idx, &col) in row.col_indices().iter().enumerate() {
                let val = Signed::abs(&row.values()[idx]);
                if col == i {
                    diag_val = val;
                } else {
                    off_diag_sum += val;
                }
            }

            if diag_val <= off_diag_sum {
                return false;
            }
        }

        true
    }

    fn scale_rows(&mut self, scaling: &DVector<T>) -> Result<()> {
        if scaling.len() != self.nrows() {
            return Err(Error::InvalidConfiguration(
                "Scaling vector size mismatch".to_string(),
            ));
        }

        // This requires row-wise access which CSR provides
        // But modification is complex, would need to rebuild
        Err(Error::InvalidConfiguration(
            "Row scaling not directly supported for CSR format".to_string(),
        ))
    }

    fn scale_columns(&mut self, scaling: &DVector<T>) -> Result<()> {
        if scaling.len() != self.ncols() {
            return Err(Error::InvalidConfiguration(
                "Scaling vector size mismatch".to_string(),
            ));
        }

        // Column scaling in CSR format - need to iterate carefully
        // Get the column indices first, then update values
        let n_entries = self.nnz();
        for idx in 0..n_entries {
            let col = self.col_indices()[idx];
            self.values_mut()[idx] *= scaling[col];
        }

        Ok(())
    }
}
