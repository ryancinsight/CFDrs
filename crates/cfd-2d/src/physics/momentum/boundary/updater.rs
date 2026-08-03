use crate::scalar::Cfd2dScalar;
use cfd_math::sparse::SparseMatrixBuilder;

/// Abstract interface for sparse matrix topological assembly and zero-allocation value updates.
///
/// Permits generalized solvers to reuse matrix buffers regardless of whether they are initially
/// constructing the CSR topology (via `SparseMatrixBuilder`) or performing in-place numerical
/// updates on an existing structure (via `SparseMatrix`).
pub trait MatrixUpdater<T> {
    /// Inserts or aggregates `val` into the block matrix at location `(row, col)`.
    ///
    /// Implementations must uphold atomic accumulation commutativity for thread-parallel execution.
    fn add_entry(&mut self, row: usize, col: usize, val: T) -> cfd_core::error::Result<()>;
}

impl<T: Cfd2dScalar + Copy> MatrixUpdater<T> for SparseMatrixBuilder<T> {
    fn add_entry(&mut self, row: usize, col: usize, val: T) -> cfd_core::error::Result<()> {
        self.add_entry(row, col, val)
    }
}

impl<T: Cfd2dScalar + Copy> MatrixUpdater<T> for cfd_math::sparse::SparseMatrix<T> {
    fn add_entry(&mut self, row: usize, col: usize, val: T) -> cfd_core::error::Result<()> {
        let row_ptr = self.row_ptr().to_vec();
        let col_indices = self.col_indices().to_vec();
        let start = row_ptr[row];
        let end = row_ptr[row + 1];
        if let Ok(idx) = col_indices[start..end].binary_search(&col) {
            self.values_mut()[start + idx] += val;
        }
        Ok(())
    }
}
