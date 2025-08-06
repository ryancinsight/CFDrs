//! Sparse matrix utilities.

use nalgebra::RealField;

/// Sparse matrix wrapper
pub type SparseMatrix<T> = nalgebra_sparse::CsrMatrix<T>;

/// Sparse matrix builder
pub struct SparseMatrixBuilder<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> SparseMatrixBuilder<T> {
    /// Create a new builder
    pub fn new(_rows: usize, _cols: usize) -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }

    /// Build the sparse matrix
    pub fn build(self) -> SparseMatrix<T> {
        nalgebra_sparse::CsrMatrix::zeros(0, 0)
    }
}